#pragma once

#include <execinfo.h>
#include <execution/sql/memory_pool.h>
#include <loggers/index_logger.h>

#include <functional>
#include <queue>
#include <shared_mutex>
#include <stack>
#include <utility>
#include <vector>

#include "common/spin_latch.h"

namespace terrier::storage::index {
/**
 * These macros are private to this file - see the #undef at the bottom
 * Convenience macro for data structure invariants. Checks whether or not
 * the condition x evaluates to true, and if it does not, it returns
 */
#define CHECK(x)                                                         \
  do {                                                                   \
    if (!(x)) {                                                          \
      INDEX_LOG_ERROR("Failed check (%s) on line: %d.\n", #x, __LINE__); \
      return false;                                                      \
    }                                                                    \
  } while (0)

/**
 * Convenience version of CHECK that checks whether x <= y
 */
#define CHECK_LE(x, y) CHECK((x) <= (y))

/**
 * Convenience version of CHECK that checks whether x < y
 */
#define CHECK_LT(x, y) CHECK((x) < (y))

/**
 * Convenience version of CHECK that checks whether x == y
 */
#define CHECK_EQ(x, y) CHECK((x) == (y))

#ifdef NDEBUG
#define DEBUG_ONLY_DEFINE(x) ((void)0)
#define DEBUG_ONLY_RUN(x) ((void)0)
#else
#define DEBUG_ONLY_DEFINE(x) x
#define DEBUG_ONLY_RUN(x) \
  do {                    \
    x                     \
  } while (0)
#endif

// Only define if you wish IsBplusTree to be run before every insert/scan
// #define DEEP_DEBUG

/**
 * Number of children an interior node is allowed to have, or
 * equivalently the number of keys a leaf node is allowed to store
 */
constexpr uint32_t NUM_CHILDREN = 256;

/**
 * Number of keys stored in an overflow node.
 */
constexpr uint32_t OVERFLOW_SIZE = 512;

/**
 * Length of time to try to increment before dying
 */
constexpr std::chrono::duration TIMEOUT = std::chrono::nanoseconds(100);

/**
 * Minimum number of children an interior node is allowed to have
 * (equivalently, the minimum number of keys a leaf node is allowed to store)
 */
constexpr uint32_t MIN_CHILDREN = NUM_CHILDREN / 2 + NUM_CHILDREN % 2;

#ifndef NDEBUG
/**
 * A lock for use in debug mode to make sure lock invariants are held.
 */
class DebugLock : public std::shared_timed_mutex {
 public:
  /**
   * The writer thread holding the lock
   */
  std::atomic<std::thread::id> holder_;
  /**
   * The number of readers holding this lock
   */
  std::atomic<uint32_t> shared_holders_;

  /**
   * A stored backtrace from the writer who last locked this
   */
  void *trace_[10];

  /**
   * The length of the backtrace of the writer who last locked this.
   */
  std::atomic<size_t> size_;
  void lock_shared() {  // NOLINT for name compatibility
    std::shared_timed_mutex::lock_shared();
    TERRIER_ASSERT(holder_ == std::thread::id(), "No thread should be a writer");
    shared_holders_++;
  }

  template <typename R, typename P>
  bool try_lock_shared_for(const std::chrono::duration<R, P> &duration) {  // NOLINT for name compatibility
    bool result = std::shared_timed_mutex::try_lock_shared_for(duration);
    if (result) {
      TERRIER_ASSERT(holder_ == std::thread::id(), "No thread should be a writer");
      shared_holders_++;
    }
    return result;
  }

  void unlock_shared() {  // NOLINT for name compatibility
    TERRIER_ASSERT(shared_holders_ > 0, "Must have at least one reader");
    shared_holders_--;
    std::shared_timed_mutex::unlock_shared();
  }

  void lock() {  // NOLINT for lock name compatibility
    if (holder_ == std::this_thread::get_id()) {
      StackTrace();
      TERRIER_ASSERT(false, "Should not lock twice!");
    }
    std::shared_timed_mutex::lock();
    TERRIER_ASSERT(shared_holders_ == 0, "Should be no readers now");
    holder_ = std::this_thread::get_id();
    size_ = backtrace(trace_, 10);
  }

  void unlock() {  // NOLINT for lock name compatibility
    TERRIER_ASSERT(holder_ == std::this_thread::get_id(), "Should only unlock from the current thread");
    holder_ = std::thread::id();
    std::shared_timed_mutex::unlock();
  }

  /**
   * Prints a stack trace from the last writer locking
   */
  void StackTrace() { backtrace_symbols_fd(trace_, size_, STDERR_FILENO); }
};
#endif  // NDEBUG

template <typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType>,
          typename KeyEqualityChecker = std::equal_to<KeyType>, typename KeyHashFunc = std::hash<KeyType>,
          typename ValueEqualityChecker = std::equal_to<ValueType>>
class BPlusTree {
 private:
  static execution::sql::MemoryPool mem_pool;

  /**
   * This inner class defines what an overflow node looks like. An overflow node
   * is defined to only contain duplicate keys of the LeafNode it is attached to
   */
  class OverflowNode {
   private:
    uint32_t filled_keys_;
    KeyType keys_[OVERFLOW_SIZE];
    ValueType values_[OVERFLOW_SIZE];
    OverflowNode *next_;

    friend class BPlusTree;

    /**
     * Creates a new overflow node from the memory pool.
     */
    static OverflowNode *CreateNew() {
      auto *node = reinterpret_cast<OverflowNode *>(mem_pool.Allocate(sizeof(OverflowNode), true));
      return node;
    }

    /**
     * Releases an overflow node back into the memory pool.
     */
    static void Delete(OverflowNode *elem) { mem_pool.Deallocate(elem, sizeof(OverflowNode)); }
  };

  /**
   * A  node represents the common portion of both a LeafNode and an InteriorNode.
   * This allows us to write code that is agnostic to the specifics of a leaf node vs an interior node,
   * but which only operates on keys and a  value.
   */
  template <typename Value>
  class GenericNode {
   protected:
    uint32_t filled_keys_ = 0;
    KeyType keys_[NUM_CHILDREN];
    Value values_[NUM_CHILDREN];
#ifndef NDEBUG
    mutable DebugLock latch_;
#else
    mutable std::shared_timed_mutex latch_;
#endif

    friend class BPlusTree;
  };

  /**
   * This inner class defines what a LeafNode is. The keys are stored in
   * sorted order with no duplicates within the keys_ array. If the
   * tree supports duplicate keys, then any duplicates will be stored
   * in the overflow node.
   */
  class LeafNode : public GenericNode<ValueType> {
   private:
    LeafNode *prev_;
    LeafNode *next_;
    OverflowNode *overflow_;

    bool Contains(KeyType k, const BPlusTree *parent) const {
      for (uint32_t i = 0; i < this->filled_keys_; i++) {
        if (parent->KeyCmpEqual(this->keys_[i], k)) {
          return true;
        }
      }
      return false;
    }

    /**
     * Checks to make sure that the given node is a leaf.
     * @param allow_duplicates If duplicates are allowed in this leaf node
     * @param is_root If this leaf node is the only leaf node in the tree
     * @param parent The BPlusTree that this leaf node is a part of
     * @return whether this node is a valid leaf
     */
    bool IsLeafNode(bool allow_duplicates, bool is_root, const BPlusTree *parent) {
      if (is_root) {
        CHECK_LE(0, this->filled_keys_);  // Root may even be empty
      } else {
        CHECK_LE(MIN_CHILDREN, this->filled_keys_);
      }
      CHECK_LE(this->filled_keys_, NUM_CHILDREN);

      CHECK(allow_duplicates || overflow_ == nullptr);

      for (uint32_t i = 1; i < this->filled_keys_; i++) {
        CHECK(parent->KeyCmpLess(this->keys_[i - 1], this->keys_[i]));
        CHECK(!parent->KeyCmpEqual(this->keys_[i - 1], this->keys_[i]));
      }

      for (OverflowNode *current = overflow_; current != nullptr; current = current->next_) {
        if (current->next_ != nullptr) {
          CHECK(current->filled_keys_ == OVERFLOW_SIZE);
        }
        for (uint32_t i = 0; i < current->filled_keys_; i++) {
          CHECK(Contains(current->keys_[i], parent));
        }
      }
      return true;
    }

    /**
     * Gets the amount of heap space used by this node.
     */
    size_t GetHeapUsage() {
      size_t usage = sizeof(LeafNode);

      for (OverflowNode *current = overflow_; current != nullptr; current = current->next_) {
        usage += sizeof(OverflowNode);
      }
      return usage;
    }

    friend class BPlusTree;

    /**
     * Creates a new leaf node from the memory pool.
     */
    static LeafNode *CreateNew() {
      auto *node = reinterpret_cast<LeafNode *>(mem_pool.Allocate(sizeof(class LeafNode), true));
      // Initialize the shared mutex
      new (&node->latch_) std::shared_timed_mutex();
      return node;
    }

    /**
     * Releases a leaf node to the memory pool.
     */
    static void Delete(LeafNode *elem) { mem_pool.Deallocate(elem, sizeof(class LeafNode)); }
  };

  class InteriorNode;

  /**
   * The possible children of an interior node.
   */
  union Child {
    InteriorNode *as_interior_;
    LeafNode *as_leaf_;
  };

  /**
   * This inner class defines what an InteriorNode is. An interior node is made up
   * of NUM_CHILDREN children, and NUM_CHILDREN - 1 guideposts between those children.
   * The children may either be other interior nodes, or leaf nodes (as dictated by the leaf_children_
   * parameter)
   *
   * NB: We allocate one more guide post than necessary in order to allow InteriorNode
   * and LeafNode to have the same layout in memory for the first 3 fields (e.g. GenericNode)
   */
  class InteriorNode : public GenericNode<Child> {
   private:
    bool leaf_children_;

    /**
     * Gets a reference to the i^th child assuming that the children of this node
     * are leaves
     */
    LeafNode *&Leaf(uint32_t i) {
      TERRIER_ASSERT(leaf_children_, "Leaf must be called only on an interior node with leaf children");
      return this->values_[i].as_leaf_;
    }

    /**
     * Gets a reference to the i^th child assuming that the children of this node
     * are interior nodes
     */
    InteriorNode *&Interior(uint32_t i) {
      TERRIER_ASSERT(!leaf_children_, "Interior must be called only on an interior node with no leaf children");
      return this->values_[i].as_interior_;
    }

    /**
     * Returns whether or not this node is an interior node
     *
     * @param allow_duplicates Does the BPlusTree allow duplicates?
     * @param prev the left sibling of this interior node, or {@code nullptr} if no left sibling exists
     * @param next the right sibling of this interior node, or {@code nullptr} if no right sibling exists
     */
    bool IsInteriorNode(bool allow_duplicates, InteriorNode *prev, InteriorNode *next, bool is_root,
                        const BPlusTree *parent) {
      // Check bounds on number of filled guide posts
      if (is_root) {
        CHECK_LT(1, this->filled_keys_);
      } else {
        CHECK_LE(MIN_CHILDREN, this->filled_keys_);
      }
      CHECK_LE(this->filled_keys_, NUM_CHILDREN);

      // Check to make sure that every child is a valid node
      for (uint32_t i = 0; i < this->filled_keys_; i++) {
        if (leaf_children_) {
          CHECK(Leaf(i) != nullptr);
          CHECK(Leaf(i)->IsLeafNode(allow_duplicates, false, parent));

          CHECK(i == 0 || (Leaf(i)->prev_ == Leaf(i - 1)));
          CHECK(i == 0 || (Leaf(i - 1)->next_ == Leaf(i)));
        } else {
          CHECK(Interior(i) != nullptr);
          InteriorNode *new_prev = nullptr;
          if (i == 0 && prev != nullptr) {
            new_prev = prev->Interior(prev->filled_keys_ - 1);
          } else if (i != 0) {
            new_prev = Interior(i - 1);
          }

          InteriorNode *new_next = nullptr;
          if (i == this->filled_keys_ - 1 && next != nullptr) {
            new_next = next->Interior(0);
          } else if (i != this->filled_keys_ - 1) {
            new_next = Interior(i + 1);
          }
          CHECK(Interior(i)->IsInteriorNode(true, new_prev, new_next, false, parent));
        }
      }

      // Check to make sure 0'th key is never touched - only if in debug mode
      /*
      DEBUG_ONLY_RUN(char *empty_key = reinterpret_cast<char *>(&this->keys_[0]);
                     for (uint32_t i = 0; i < sizeof(KeyType); i++) { CHECK_EQ(empty_key[i], '\0'); });
                     */

      // Make sure guide posts are in sorted order with no dupes
      for (uint32_t i = 2; i < this->filled_keys_; i++) {
        CHECK(parent->KeyCmpLess(this->keys_[i - 1], this->keys_[i]));
      }

      // Check to make sure that each child has keys that are in the correct range
      for (uint32_t i = 1; i < this->filled_keys_; i++) {
        if (leaf_children_) {
          auto leaf = Leaf(i);
          auto prev_leaf = Leaf(i - 1);
          CHECK(parent->KeyCmpLessEqual(this->keys_[i], leaf->keys_[0]));
          CHECK(parent->KeyCmpLessEqual(prev_leaf->keys_[prev_leaf->filled_keys_ - 1], this->keys_[i]));
        } else {
          auto interior = Interior(i);
          while (!interior->leaf_children_) {
            interior = interior->Interior(0);
          }
          auto prev_interior = Interior(i - 1);
          CHECK(parent->KeyCmpLessEqual(this->keys_[i], interior->Leaf(0)->keys_[0]));
          CHECK(parent->KeyCmpLessEqual(prev_interior->keys_[prev_interior->filled_keys_ - 1], this->keys_[i]));
        }
      }

      // Check to make sure that children on the edges (e.g. at indices 0 and filled_keys_)
      // have keys in the correct ranges
      if (leaf_children_) {
        CHECK(prev == nullptr || prev->leaf_children_);
        CHECK((prev == nullptr && Leaf(0)->prev_ == nullptr) ||
              (prev != nullptr && Leaf(0)->prev_ == prev->Leaf(prev->filled_keys_ - 1)));
        CHECK(prev == nullptr || parent->KeyCmpLessEqual(prev->keys_[prev->filled_keys_ - 1], Leaf(0)->keys_[0]));

        auto last_leaf = Leaf(this->filled_keys_ - 1);
        CHECK(next == nullptr || next->leaf_children_);
        CHECK((next == nullptr && last_leaf->next_ == nullptr) ||
              (next != nullptr && last_leaf->next_ == next->Leaf(0)));

        CHECK(next == nullptr ||
              parent->KeyCmpLessEqual(last_leaf->keys_[last_leaf->filled_keys_ - 1], next->Leaf(0)->keys_[0]));
      } else {
        CHECK(prev == nullptr || parent->KeyCmpLessEqual(prev->keys_[prev->filled_keys_ - 1], Interior(0)->keys_[1]));
        auto last_interior = Interior(this->filled_keys_ - 1);
        CHECK(next == nullptr ||
              parent->KeyCmpLessEqual(last_interior->keys_[last_interior->filled_keys_ - 1], next->keys_[1]));
      }

      return true;
    }

    /**
     * Gets the heap usage of this node, in bytes.
     */
    size_t GetHeapUsage() {
      size_t usage = sizeof(InteriorNode);

      for (uint32_t i = 0; i < this->filled_keys_; i++) {
        if (leaf_children_) {
          usage += Leaf(i)->GetHeapUsage();
        } else {
          usage += Interior(i)->GetHeapUsage();
        }
      }

      return usage;
    }

    friend class BPlusTree;

    /**
     * Creates a new InteriorNode from the memory pool.
     */
    static InteriorNode *CreateNew() {
      auto *node = reinterpret_cast<InteriorNode *>(mem_pool.Allocate(sizeof(class InteriorNode), true));
      // Initialize the shared mutex
      new (&node->latch_) std::shared_timed_mutex();
      return node;
    }

    /**
     * Releases an InteriorNode back into the memory pool.
     */
    static void Delete(InteriorNode *elem) { mem_pool.Deallocate(elem, sizeof(class InteriorNode)); }
  };

  InteriorNode *root_;
  uint32_t depth_;
  bool allow_duplicates_;
  KeyComparator key_cmp_obj_;
  KeyEqualityChecker key_eq_obj_;
  KeyHashFunc key_hash_obj_;
  ValueEqualityChecker value_eq_obj_;

  /**
   * Inserts the key/value pair into a node that still has room
   */
  template <typename Value>
  void InsertIntoNode(GenericNode<Value> *node, uint32_t insertionSpot, KeyType k, Value v) {
    TERRIER_ASSERT(node->filled_keys_ < NUM_CHILDREN, "There must still be room in the node!");

    for (uint32_t j = node->filled_keys_; j > insertionSpot; j--) {
      node->keys_[j] = node->keys_[j - 1];
      node->values_[j] = node->values_[j - 1];
    }

    node->keys_[insertionSpot] = k;
    node->values_[insertionSpot] = v;
    node->filled_keys_++;
  }

  /**
   * Splits a node to allow for the insertion of a key/value pair. This method will split the node {@code from}
   * into both {@code from} and {@code to}, where {@code from} will be the left sibling of {@code to}
   *
   * @param from The node to split in half
   * @param to The node that we will insert into
   * @param insertion_spot The location in {@code from} that the key/value pair is supposed to be inserted
   * @param k the key to insert
   * @param v the value to insert
   * @return A guide post that splits from and to such that every element in from is less than or equal to every element
   *         in to.
   */
  template <typename Value>
  KeyType SplitNode(GenericNode<Value> *from, GenericNode<Value> *to, uint32_t insertion_spot, const KeyType &k,
                    const Value &v) {
    TERRIER_ASSERT(from != nullptr && to != nullptr, "No nullptrs allowed here!");
    TERRIER_ASSERT(from->filled_keys_ == NUM_CHILDREN, "Split node should only be called on a full node!");
    // NB: This assert cannot properly check for the node being completely empty, since an interior node is
    //     full whenever the number of "filled keys" is 1 and a leaf node when the number of filled keys is 0
    TERRIER_ASSERT(to->filled_keys_ == 0 || to->filled_keys_ == 1,
                   "Split node should be called to split into an empty node!");

    const KeyType *guide_post_ptr = nullptr;
    if (insertion_spot <= NUM_CHILDREN / 2) {
      guide_post_ptr = &(from->keys_[NUM_CHILDREN / 2]);
    } else if (insertion_spot == NUM_CHILDREN / 2 + 1) {
      guide_post_ptr = &k;
    } else {
      guide_post_ptr = &(from->keys_[NUM_CHILDREN / 2 + 1]);
    }
    TERRIER_ASSERT(guide_post_ptr != nullptr, "Guide post should not be null!");
    const KeyType guide_post = *guide_post_ptr;

    DEBUG_ONLY_DEFINE(bool is_interior = to->filled_keys_ == 1);
    to->filled_keys_ = 0;
    if (insertion_spot <= NUM_CHILDREN / 2) {
      for (uint32_t j = NUM_CHILDREN / 2; j < NUM_CHILDREN; ++j) {
        to->keys_[to->filled_keys_] = from->keys_[j];
        to->values_[to->filled_keys_] = from->values_[j];
        ++to->filled_keys_;
      }

      from->filled_keys_ = NUM_CHILDREN / 2;
      InsertIntoNode(from, insertion_spot, k, v);
    } else {
      for (uint32_t j = NUM_CHILDREN / 2 + 1; j < insertion_spot; ++j) {
        to->keys_[to->filled_keys_] = from->keys_[j];
        to->values_[to->filled_keys_] = from->values_[j];
        ++to->filled_keys_;
      }

      to->keys_[to->filled_keys_] = k;
      to->values_[to->filled_keys_] = v;
      ++to->filled_keys_;

      for (uint32_t j = insertion_spot; j < NUM_CHILDREN; ++j) {
        to->keys_[to->filled_keys_] = from->keys_[j];
        to->values_[to->filled_keys_] = from->values_[j];
        ++to->filled_keys_;
      }

      from->filled_keys_ = NUM_CHILDREN / 2 + 1;
    }

    DEBUG_ONLY_RUN(if (is_interior) { memset(&to->keys_[0], 0, sizeof(KeyType)); });

    return guide_post;
  }

  /**
   * Finds the first index i in node where node->keys_[i] >= k,
   * or node->filled_keys_ if such an index does not exist
   */
  uint32_t FindKey(const LeafNode *node, KeyType k) const {
    uint32_t lo = 0;
    uint32_t hi = node->filled_keys_;
    while (lo < hi) {
      uint32_t mid = (lo + hi) / 2;
      if (KeyCmpLessEqual(k, node->keys_[mid])) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  /**
   * Finds the first index i in node where node->keys_[i] > k,
   * or node->filled_keys_ if such a node does not exist.
   */
  uint32_t FindKey(const InteriorNode *node, KeyType k) const {
    uint32_t lo = 1;
    uint32_t hi = node->filled_keys_;
    while (lo < hi) {
      uint32_t mid = (lo + hi) / 2;
      if (KeyCmpLess(k, node->keys_[mid])) {
        hi = mid;
      } else {
        lo = mid + 1;
      }
    }
    return lo;
  }

  /**
   * Traverses the path from the root to the leaf that should contain k (if k exists in this tree).
   * Keeps track of the interior nodes that may need to be modified in case of a split in the leaf in
   * the `potential_changes` vector.
   *
   * @return the leaf that should contain k
   */
  LeafNode *TraverseTrack(InteriorNode *root, KeyType k, std::vector<InteriorNode *> *potential_changes) const {
    InteriorNode *current = root;
    /**
    uint32_t max_depth = depth_;
    uint32_t cur_depth = 1;*/
    LeafNode *leaf = nullptr;
    while (leaf == nullptr) {
      // If we know we do not need to split (e.g. there are still open guide posts,
      // then we don't need to keep track of anything above us
      if (potential_changes != nullptr && current->filled_keys_ < NUM_CHILDREN) {
        for (const auto &interior : *potential_changes) {
          interior->latch_.unlock();
        }
        potential_changes->erase(potential_changes->begin(), potential_changes->end());
      }

      // However, the level below us may still split, so we may still change. Keep track of that
      if (potential_changes != nullptr) {
        potential_changes->push_back(current);
      }

      uint32_t i = FindKey(current, k) - 1;

      // If we've reached the leaf level, break out!
      if (current->leaf_children_) {
        leaf = current->Leaf(i);
        leaf->latch_.lock();
        if (potential_changes == nullptr) {
          current->latch_.unlock_shared();
        }
        TERRIER_ASSERT(leaf != nullptr, "Leaves should not be null!!");
      } else {
        InteriorNode *prev = current;
        current = current->Interior(i);
        if (potential_changes == nullptr) {
          current->latch_.lock_shared();
          prev->latch_.unlock_shared();
        } else {
          current->latch_.lock();
        }
      }
      // cur_depth++;
    }

    // If the leaf definitely has enough room, then we don't need to split anything!
    if (potential_changes != nullptr && leaf->filled_keys_ < NUM_CHILDREN) {
      for (const auto &interior : *potential_changes) {
        interior->latch_.unlock();
      }
      potential_changes->erase(potential_changes->begin(), potential_changes->end());
    }

    return leaf;
  }

  /**
   * Splits a full leaf into two leaves and inserts k and v into the appropriate slot
   */
  KeyType SplitLeaf(LeafNode *from, LeafNode *to, uint32_t insertion_spot, KeyType k, ValueType v) {
    KeyType guide_post = SplitNode(reinterpret_cast<GenericNode<ValueType> *>(from),
                                   reinterpret_cast<GenericNode<ValueType> *>(to), insertion_spot, k, v);

    // These two variables represent where in from's overflow nodes we are reading from
    OverflowNode *read_overflow = from->overflow_;
    uint32_t read_id = 0;

    // These two variables represent where in from's overflow nodes we are writing to
    OverflowNode *write_overflow = from->overflow_;
    uint32_t write_id = 0;

    // These two variables represent where in to's overflow nodes we are writing to
    OverflowNode *to_overflow = nullptr;
    uint32_t to_id = 0;

    // Unlike in from, we might have to allocate new overflow leaves. This is where we should write the
    // pointer in order to ensure that the link is created - we do it this way to support lazy allocation
    OverflowNode **to_overflow_alloc = &(to->overflow_);

    // Skip over any overflow elements that don't need to be moved at all
    while (read_overflow != nullptr && KeyCmpLess(read_overflow->keys_[read_id], guide_post)) {
      read_id++;
      write_id++;
      while (read_id == OVERFLOW_SIZE || (read_overflow != nullptr && read_id >= read_overflow->filled_keys_)) {
        read_overflow = read_overflow->next_;
        read_id = 0;
        write_overflow = write_overflow->next_;
        write_id = 0;
      }
    }

    if (read_overflow != nullptr && read_id >= read_overflow->filled_keys_) {
      read_overflow = read_overflow->next_;
      read_id = 0;
    }

    // For every other overflow element, move it either from from to from (in another slot) or from from to
    // new from (in a new slot) depending on where it should live
    while (read_overflow != nullptr) {
      if (KeyCmpLess(read_overflow->keys_[read_id], guide_post)) {
        if (write_id == OVERFLOW_SIZE) {
          write_overflow = write_overflow->next_;
          write_id = 0;
        }
        // Move from from to from
        write_overflow->keys_[write_id] = read_overflow->keys_[read_id];
        write_overflow->values_[write_id] = read_overflow->values_[read_id];

        // Advance the write index
        write_id++;
      } else {
        // Lazy allocation still means we have to allocate eventually, and now we know we need allocate
        if (to_overflow == nullptr) {
          to_overflow = OverflowNode::CreateNew();
          *to_overflow_alloc = to_overflow;
          to_overflow_alloc = &(to_overflow->next_);
        }

        // Move from from to to.
        to_overflow->keys_[to_id] = read_overflow->keys_[read_id];
        to_overflow->values_[to_id] = read_overflow->values_[read_id];

        // Advance the write index
        to_id++;
        if (to_id == OVERFLOW_SIZE) {
          to_overflow->filled_keys_ = OVERFLOW_SIZE;
          to_overflow = nullptr;
          to_id = 0;
        }
      }

      // Advance the read index
      read_id++;
      while (read_id == OVERFLOW_SIZE || (read_overflow != nullptr && read_id >= read_overflow->filled_keys_)) {
        read_overflow = read_overflow->next_;
        read_id = 0;
      }
    }

    // If we wrote any values to the to, we now might need to
    // free some overflow nodes!
    if (write_overflow != nullptr) {
      write_overflow->filled_keys_ = write_id;
      OverflowNode *temp = write_overflow;
      write_overflow = write_overflow->next_;
      temp->next_ = nullptr;
      while (write_overflow != nullptr) {
        temp = write_overflow->next_;
        OverflowNode::Delete(write_overflow);
        write_overflow = temp;
      }
    }

    // If we did not completely fill up the new from overflow nodes, mark where we filled up to
    if (to_overflow != nullptr) {
      to_overflow->filled_keys_ = to_id;
    }

    to->next_ = from->next_;
    from->next_ = to;
    to->prev_ = from;
    if (to->next_ != nullptr) {
      to->next_->prev_ = to;
    }

    TERRIER_ASSERT(from->IsLeafNode(allow_duplicates_, false, this), "Old leaf was not preserved as leaf");
    TERRIER_ASSERT(to->IsLeafNode(allow_duplicates_, false, this), "New Leaf not preserved as leaf");

    from->latch_.unlock();
    return guide_post;
  }

  /**
   * Splits a full interior node into two interior nodes and inserts k and v into the appropriate slot
   */
  KeyType SplitInterior(InteriorNode *from, InteriorNode *to, KeyType guide_post, InteriorNode *child) {
    // Find out where to insert this node
    uint32_t i = FindKey(from, guide_post);
    TERRIER_ASSERT(i >= from->filled_keys_ || !KeyCmpEqual(guide_post, from->keys_[i]),
                   "We should not have duplicated guide posts!");

    // Perform the insert
    to->filled_keys_ = 1;
    guide_post = SplitNode(reinterpret_cast<GenericNode<InteriorNode *> *>(from),
                           reinterpret_cast<GenericNode<InteriorNode *> *>(to), i, guide_post, child);

    TERRIER_ASSERT(from->values_[0].as_interior_ != nullptr, "First value should not be nullptr!");
    TERRIER_ASSERT(to->values_[0].as_interior_ != nullptr, "First value should not be nullptr!");
    return guide_post;
  }

 public:
  /**
   * Creates the BPlusTree as an empty tree, with the given comparators
   * @param allow_duplicates Whether this tree should allow duplicatees.
   */
  explicit BPlusTree(bool allow_duplicates, KeyComparator key_cmp_obj = KeyComparator{},
                     KeyEqualityChecker key_eq_obj = KeyEqualityChecker{}, KeyHashFunc key_hash_obj = KeyHashFunc{},
                     ValueEqualityChecker value_eq_obj = ValueEqualityChecker{})
      : depth_(1),
        allow_duplicates_(allow_duplicates),
        key_cmp_obj_(key_cmp_obj),
        key_eq_obj_(key_eq_obj),
        key_hash_obj_(key_hash_obj),
        value_eq_obj_(value_eq_obj) {
    root_ = InteriorNode::CreateNew();
    root_->leaf_children_ = true;
    root_->filled_keys_ = 1;
    root_->Leaf(0) = LeafNode::CreateNew();
  }

  /**
   * Checks if this B+ Tree is truly a B+ Tree
   */
  [[nodiscard]] bool IsBplusTree() const {
    if (root_->filled_keys_ == 1) {
      // Root is not a valid interior node until we do a split!
      // However, this is only allowed to happen if the depth is truly 1
      CHECK_EQ(depth_, 1);
      CHECK(root_->leaf_children_);
      CHECK(root_->Leaf(0)->IsLeafNode(allow_duplicates_, true, this));
      CHECK(root_->Leaf(0)->next_ == nullptr && root_->Leaf(0)->prev_ == nullptr);
      return true;
    }

    if (!root_->IsInteriorNode(allow_duplicates_, nullptr, nullptr, true, this)) {
      return false;
    }

    // Check to make sure that the depth of the tree is the same regardless of
    // which way you go down the tree
    std::queue<InteriorNode *> level;
    level.push(root_);
    for (uint32_t i = 1; i < depth_; i++) {
      std::queue<InteriorNode *> next_level;
      while (!level.empty()) {
        auto elem = level.front();
        level.pop();
        CHECK(!elem->leaf_children_);
        for (uint32_t j = 0; j < elem->filled_keys_; j++) {
          next_level.push(elem->Interior(j));
        }
      }
      level = std::move(next_level);
    }

    while (!level.empty()) {
      auto elem = level.front();
      level.pop();
      CHECK(elem->leaf_children_);
    }
    return true;
  }

  /*
   * KeyCmpLess() - Compare two keys for "less than" relation
   *
   * If key1 < key2 return true
   * If not return false
   */
  bool KeyCmpLess(const KeyType &key1, const KeyType &key2) const { return key_cmp_obj_(key1, key2); }

  /*
   * KeyCmpEqual() - Compare a pair of keys for equality
   *
   * This functions compares keys for equality relation
   */
  bool KeyCmpEqual(const KeyType &key1, const KeyType &key2) const { return key_eq_obj_(key1, key2); }

  /*
   * KeyCmpGreaterEqual() - Compare a pair of keys for >= relation
   *
   * It negates result of keyCmpLess()
   */
  bool KeyCmpGreaterEqual(const KeyType &key1, const KeyType &key2) const { return !KeyCmpLess(key1, key2); }

  /*
   * KeyCmpGreater() - Compare a pair of keys for > relation
   *
   * It flips input for keyCmpLess()
   */
  bool KeyCmpGreater(const KeyType &key1, const KeyType &key2) const { return KeyCmpLess(key2, key1); }

  /*
   * KeyCmpLessEqual() - Compare a pair of keys for <= relation
   */
  bool KeyCmpLessEqual(const KeyType &key1, const KeyType &key2) const { return !KeyCmpGreater(key1, key2); }

  // Forward declaration for friend classing
  class KeyIterator;

  /**
   * An iterator used for iterating through the multiple values which could exist for a single key.
   */
  class DuplicateIterator {
   private:
    friend class KeyIterator;

    const BPlusTree *tree_;
    OverflowNode *current_;
    uint32_t index_;

    DuplicateIterator(const BPlusTree *tree, OverflowNode *current, uint32_t index)
        : tree_(tree), current_(current), index_(index) {}

   public:
    /**
     * Grabs the value currently referenced by this iterator.
     */
    ValueType &operator*() const {
      TERRIER_ASSERT(current_ != nullptr, "Derefernce ofinvalid iterator");
      return current_->values_[index_];
    }

    /**
     * Moves to the next value, or the end if there are no more values.
     * @returns a reference to this iterator.
     */
    DuplicateIterator &operator++() {
      KeyType &key = current_->keys_[index_];
      do {
        index_++;
        if (index_ >= current_->filled_keys_) {
          current_ = current_->next_;
          index_ = 0;
        }
      } while (current_ != nullptr && !tree_->KeyCmpEqual(key, current_->keys_[index_]));
      return *this;
    }

    /**
     * Moves to the next value, or the end if there are no more values.
     * @returns a copy of the iterator before the operation.
     */
    DuplicateIterator operator++(int) {
      KeyIterator copy = *this;
      this->operator++();
      return copy;
    }

    /**
     * Checks if two iterators are referring to the same duplicate value.
     */
    bool operator==(const DuplicateIterator &other) const {
      return current_ == other.current_ && index_ == other.index_;
    }

    /**
     * Checks if two iterators are not referring to the same duplicate value.
     */
    bool operator!=(const DuplicateIterator &other) const { return !this->operator==(other); }
  };

  /**
   * An iterator over specific keys in the BPlusTree. Use ++ and -- to move to greater / less keys.
   *
   * Will only visit unique keys in the tree. Duplicate values can be gotten by using a duplicate iterator,
   * gotten by this one's begin() method.
   */
  class KeyIterator {
   private:
    friend class BPlusTree;

    const BPlusTree *tree_;
    LeafNode *current_;
    uint32_t index_;
    bool valid_;

    KeyIterator(const BPlusTree *tree, LeafNode *current, uint32_t index)
        : tree_(tree), current_(current), index_(index), valid_(true) {}

   public:
    /**
     * KeyIterators can be invalidated if they are advancing and encounter a situation
     * which would cause a deadlock. Use this to check for that and restart a scan.
     * @returns true if valid, false otherwise.
     */
    bool Valid() const { return valid_; }

    /**
     * Gets the primary (first-inserted) value stored for this key.
     */
    ValueType &Value() const {
      TERRIER_ASSERT(current_ != nullptr, "Value called on invalid iterator");
      return current_->values_[index_];
    }

    /**
     * Gets the primary (first-inserted) value stored for this key.
     *
     * Syntactic-sugar for STL compatability
     */
    ValueType &operator*() const {
      TERRIER_ASSERT(current_ != nullptr, "Value called on invalid iterator");
      return current_->values_[index_];
    }

    /**
     * Gets the key for this iterator.
     */
    const KeyType &Key() const {
      TERRIER_ASSERT(current_ != nullptr, "Key called on invalid iterator");
      return current_->keys_[index_];
    }

    /**
     * Releases the lock the iterator holds on its current key and immediately invalidates this iterator
     */
    void ReleaseLock() {
      if (current_ != nullptr) {
        current_->latch_.unlock_shared();
      }
      current_ = nullptr;
      index_ = 0;
    }

    /**
     * Moves to the immediately greater key, or the end of the tree if there are no more keys.
     * @returns a reference to this iterator.
     */
    KeyIterator &operator++() {
      index_++;
      if (index_ >= current_->filled_keys_) {
        auto *next = current_->next_;
        if (next) {
          bool acquired = next->latch_.try_lock_shared_for(TIMEOUT);
          if (!acquired) {
            valid_ = false;
            current_->latch_.unlock_shared();
            return *this;
          }
        }
        current_->latch_.unlock_shared();
        current_ = next;
        index_ = 0;
      }
      return *this;
    }

    /**
     * Moves to the immediately greater key, or the end of the tree if there are no more keys.
     * @returns a copy of the iterator before the operation.
     */
    KeyIterator operator++(int) {
      KeyIterator copy = *this;
      this->operator++();
      return copy;
    }

    /**
     * Moves to the immediately less-than key, or the end of the tree if there are no more keys.
     * @returns a reference to this iterator.
     */
    KeyIterator &operator--() {
      TERRIER_ASSERT(current_ != nullptr, "-- called on an invalid iterator!");

      if (index_ == 0) {
        auto *prev = current_->prev_;
        if (prev != nullptr) {
          bool acquired = prev->latch_.try_lock_shared_for(TIMEOUT);
          if (!acquired) {
            valid_ = false;
            current_->latch_.unlock_shared();
            return *this;
          }
        }
        current_->latch_.unlock_shared();
        current_ = prev;
        if (current_ != nullptr) {
          index_ = current_->filled_keys_ - 1;
        }
      } else {
        index_--;
      }
      return *this;
    }

    /**
     * Moves to the immediately less-than key, or the end of the tree if there are no more keys.
     * @returns a copy of the iterator before the operation.
     */
    KeyIterator operator--(int) {
      KeyIterator copy = *this;
      this->operator--();
      return copy;
    }

    /**
     * Checks if two iterators are referring to the same key.
     */
    bool operator==(const KeyIterator &other) const {
      return other.current_ == this->current_ && other.index_ == this->index_;
    }

    /**
     * Checks if two iterators are not referring to the same key.
     */
    bool operator!=(const KeyIterator &other) const { return !this->operator==(other); }

    /**
     * Gets the starting iterator for iterating over other stored values for this key.
     */
    DuplicateIterator begin() const {  // NOLINT for STL name compability
      KeyType &key = current_->keys_[index_];
      OverflowNode *current = current_->overflow_;
      uint32_t index = 0;
      while (current != nullptr && !tree_->KeyCmpEqual(key, current->keys_[index])) {
        index++;
        if (index >= current->filled_keys_) {
          current = current->next_;
          index = 0;
        }
      }
      return {tree_, current, index};
    }

    /**
     * Gets the ending iterator for iterating over other stored values for this key.
     */
    const DuplicateIterator end() const {  // NOLINT for STL name compability
      return {tree_, nullptr, 0};
    }
  };

  /**
   * Gets an iterator pointing to the first key/value pair in the tree
   * @return the iterator
   */
  KeyIterator begin() const {  // NOLINT for STL name compability
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "begin must be called on a valid B+ Tree");
#endif
    InteriorNode *save = root_;
    save->latch_.lock_shared();
    while (save != root_) {
      save->latch_.unlock_shared();
      save = root_;
      save->latch_.lock_shared();
    }

    InteriorNode *current = root_;
    while (!current->leaf_children_) {
      auto *next = current->Interior(0);
      next->latch_.lock_shared();
      current->latch_.unlock_shared();
      current = current->Interior(0);
    }

    auto *leaf = current->Leaf(0);
    leaf->latch_.lock_shared();
    current->latch_.unlock_shared();
    if (leaf->filled_keys_ == 0) {
      // Special case: empty tree
      leaf->latch_.unlock_shared();
      return {this, nullptr, 0};
    }
    return {this, leaf, 0};
  }

  /**
   * Gets an iterator pointing to the first key/value pair in the tree, whose key is >= the key passed in.
   * If the passed in key is larger than all keys in the tree, returns the end iterator.
   * @param key the Lower bound (inclusive) for the
   * @return the iterator
   */
  KeyIterator BeginGreaterEqual(const KeyType &key) const {
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "begin must be called on a valid B+ Tree");
#endif

    InteriorNode *save = root_;
    save->latch_.lock_shared();
    while (save != root_) {
      save->latch_.unlock_shared();
      save = root_;
      save->latch_.lock_shared();
    }

    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;

    // Do a search for the right-most leaf with keys < this key
    while (leaf == nullptr) {
      uint32_t child = FindKey(current, key) - 1;
      if (current->leaf_children_) {
        leaf = current->Leaf(child);
        leaf->latch_.lock_shared();
        current->latch_.unlock_shared();
        TERRIER_ASSERT(leaf != nullptr, "Leaf should be reached");
      } else {
        auto *next = current->Interior(child);
        next->latch_.lock_shared();
        current->latch_.unlock_shared();
        current = next;
      }
    }

    // Find first key in this node >= this one
    uint32_t index = FindKey(leaf, key);
    if (index < leaf->filled_keys_) {
      return {this, leaf, index};
    }

    // Key exists in the next node over
    LeafNode *next = leaf->next_;
    if (next) {
      bool acquired = next->latch_.try_lock_shared_for(TIMEOUT);
      if (!acquired) {
        KeyIterator result = {this, next, 0};
        result.valid_ = false;
        leaf->latch_.unlock_shared();
        return result;
      }
    }
    leaf->latch_.unlock_shared();
    return {this, next, 0};
  }

  /**
   * Gets an iterator pointing to the right-most value in the tree <= to the passed-in
   * key. Can be the end iterator if all keys are > than the key.
   * @param key The key to compare.
   * @returns The iterator, with a lock on the leaf (if there is one)
   */
  KeyIterator BeginLessEqual(const KeyType &key) const {
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "begin must be called on a valid B+ Tree");
#endif

    InteriorNode *save = root_;
    save->latch_.lock_shared();
    while (save != root_) {
      save->latch_.unlock_shared();
      save = root_;
      save->latch_.lock_shared();
    }

    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;

    // Do a search for the right-most leaf with keys < this key
    while (leaf == nullptr) {
      uint32_t child = FindKey(current, key) - 1;
      if (current->leaf_children_) {
        leaf = current->Leaf(child);
        leaf->latch_.lock_shared();
        current->latch_.unlock_shared();
        TERRIER_ASSERT(leaf != nullptr, "Leaf should be reached");
      } else {
        auto *next = current->Interior(child);
        next->latch_.lock_shared();
        current->latch_.unlock_shared();
        current = next;
      }
    }

    // Find last key in this node <= this one
    uint32_t index = FindKey(leaf, key);
    if (index == leaf->filled_keys_ || !KeyCmpEqual(key, leaf->keys_[index])) {
      // index should be to the left
      if (index == 0) {
        // Key exists in the leaf to the left
        KeyIterator ret(this, leaf, index);
        --ret;
        return ret;
      }
      return {this, leaf, index - 1};
    }

    return {this, leaf, index};
  }

  /**
   * Gets the end iterator.
   */
  const KeyIterator end() const {  // NOLINT for STL name compability
    return {this, nullptr, 0};
  }

  /**
   * Inserts a key/value pair into the bplustree.
   * @param k The key to insert at
   * @param v The value to insert
   * @param allow_duplicates Whether to allow duplicates for this insertion
   * @param predicate If there are no duplicates, this function can allow you to optionally
   * override a value if it returns true for that value.
   * @returns true if the insert was successful, false otherwise.
   */
  bool Insert(
      KeyType k, ValueType v, bool allow_duplicates = true,
      const std::function<bool(const ValueType &)> &predicate = [](ValueType v) { return false; }) {
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "Insert must be called on a valid B+ Tree");
#endif

    // Optimistically grab read latches on the way down
    InteriorNode *save = root_;
    save->latch_.lock_shared();
    while (save != root_) {
      save->latch_.unlock_shared();
      save = root_;
      save->latch_.lock_shared();
    }
    this->allow_duplicates_ = allow_duplicates;

    LeafNode *leaf = TraverseTrack(root_, k, nullptr);

    uint32_t i = FindKey(leaf, k);

    std::vector<InteriorNode *> potential_changes;  // Mark the interior nodes that may be split
    if ((i == leaf->filled_keys_ || !KeyCmpEqual(k, leaf->keys_[i])) && leaf->filled_keys_ == NUM_CHILDREN) {
      leaf->latch_.unlock();
      // We will need to split. guess optimism was not good enough :(
      potential_changes.reserve(depth_);

      save = root_;
      save->latch_.lock();
      while (save != root_) {
        save->latch_.unlock();
        save = root_;
        save->latch_.lock();
      }

      leaf = TraverseTrack(root_, k, &potential_changes);

      i = FindKey(leaf, k);
    }

    // If duplicates are not allowed, do not insert duplicates!
    if (!allow_duplicates && i != leaf->filled_keys_ && KeyCmpEqual(k, leaf->keys_[i])) {
      for (const auto &interior : potential_changes) {
        interior->latch_.unlock();
      }

      if (!predicate(leaf->values_[i])) {
        leaf->values_[i] = v;
        leaf->latch_.unlock();
        return true;
      }

      leaf->latch_.unlock();
      return false;
    }

    // If duplicates _are_ allowed, then insert into the overflow node, allocating a new one as necessary
    if (i != leaf->filled_keys_ && KeyCmpEqual(k, leaf->keys_[i])) {
      for (const auto &interior : potential_changes) {
        interior->latch_.unlock();
      }
      OverflowNode *current_overflow = leaf->overflow_;

      // In case at any point current_overflow is null, keep track of the location where the previous block
      // will need to point to its new "next overflow block"
      OverflowNode **prev_overflow = &(leaf->overflow_);
      while (current_overflow != nullptr && current_overflow->filled_keys_ >= OVERFLOW_SIZE) {
        prev_overflow = &(current_overflow->next_);
        current_overflow = current_overflow->next_;
      }
      if (current_overflow == nullptr) {
        current_overflow = OverflowNode::CreateNew();
        *prev_overflow = current_overflow;
      }

      current_overflow->keys_[current_overflow->filled_keys_] = k;
      current_overflow->values_[current_overflow->filled_keys_] = v;
      current_overflow->filled_keys_++;
      leaf->latch_.unlock();
      return true;
    }

    // We need to split!
    if (leaf->filled_keys_ == NUM_CHILDREN) {
      auto *new_leaf = LeafNode::CreateNew();
      KeyType guide_post = SplitLeaf(leaf, new_leaf, i, k, v);

      TERRIER_ASSERT(!potential_changes.empty(), "Potential changes should not be empty!!");
      TERRIER_ASSERT(
          potential_changes.front()->filled_keys_ < NUM_CHILDREN || potential_changes.front() == root_,
          "Potential changes should contain a front which has room for a new child, or the front should be root_");

      // to_insert represents the current node to insert into the parent.
      Child to_insert;
      to_insert.as_leaf_ = new_leaf;

      // inner represents the parent (where we insert)
      auto inner = --potential_changes.end();
      bool leaf_children = true;
      for (; inner != potential_changes.begin(); --inner) {
        auto *new_inner = InteriorNode::CreateNew();
        new_inner->leaf_children_ = leaf_children;
        leaf_children = false;
        guide_post = SplitInterior(*inner, new_inner, guide_post, to_insert.as_interior_);
        (*inner)->latch_.unlock();
        to_insert.as_interior_ = new_inner;
      }

      // Can we insert here?
      if ((*inner)->filled_keys_ < NUM_CHILDREN) {
        // Yes! Just insert!
        i = FindKey(*inner, guide_post);
        TERRIER_ASSERT(i >= (*inner)->filled_keys_ || !KeyCmpEqual(guide_post, (*inner)->keys_[i]),
                       "We should not have duplicated guide posts!");

        InsertIntoNode(*inner, i, guide_post, to_insert);
        (*inner)->latch_.unlock();
      } else {
        // No - we've reached the root and must split
        TERRIER_ASSERT(*inner == root_, "Top level inner potential change should not be full unless it is the root");

        // Split the root
        auto *new_inner = InteriorNode::CreateNew();
        new_inner->leaf_children_ = leaf_children;
        guide_post = SplitInterior(*inner, new_inner, guide_post, to_insert.as_interior_);

        auto *new_left = InteriorNode::CreateNew();
        for (uint32_t j = 0; j < root_->filled_keys_; j++) {
          new_left->keys_[j] = root_->keys_[j];
          new_left->values_[j] = root_->values_[j];
        }
        new_left->filled_keys_ = root_->filled_keys_;
        new_left->leaf_children_ = root_->leaf_children_;

        root_->leaf_children_ = false;
        root_->keys_[1] = guide_post;
        root_->filled_keys_ = 2;
        root_->Interior(0) = new_left;
        root_->Interior(1) = new_inner;
        depth_++;
        (*inner)->latch_.unlock();
      }
    } else {
      TERRIER_ASSERT(potential_changes.empty(), "Should not have any potential changes if we did not need to split!");
      InsertIntoNode(leaf, i, k, v);
      leaf->latch_.unlock();
    }

#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "Insert must return a valid B+ Tree");
#endif
    return true;
  }

  /**
   * Traverses tree to correct key, keeping track of siblings and indices needed to get to child or value.
   * @param root The root of the tree to traverse.
   * @param k The key to look for
   * @param potential_changes A vector which will be filled with potential interior nodes which need changing.
   * @param indices A vector which will be filled with various indices we take as we go down the tree.
   * @param left_siblings A vector which will be filled with the current left sibling of the interior node at
   * the equivalent node in potential_changes. Could be nullptr if at the left edge.
   * @param right_siblings A vector which will be filled with the current right sibling of the interior node at
   * the equivalent node in potential_changes. Could be nullptr if at the right edge.
   * @returns The leaf node at the bottom which we find.
   */
  LeafNode *TraverseTrackWithSiblings(InteriorNode *root, KeyType k, std::vector<InteriorNode *> *potential_changes,
                                      std::vector<uint32_t> *indices, std::vector<InteriorNode *> *left_siblings,
                                      std::vector<InteriorNode *> *right_siblings) {
    potential_changes->reserve(depth_);
    indices->reserve(depth_);
    left_siblings->reserve(depth_);
    right_siblings->reserve(depth_);

    LeafNode *leaf = nullptr;
    InteriorNode *current = root;
    InteriorNode *left = nullptr;
    InteriorNode *right = nullptr;
    uint32_t i;

    while (leaf == nullptr) {
      // This is the case where we know the current level will not need to merge. We don't need to keep track of the
      // parents anymore. We preserve the last node since future levels could potentially need it.
      if (current->filled_keys_ > MIN_CHILDREN) {
        for (const auto &interior : *potential_changes) {
          interior->latch_.unlock();
        }
        for (const auto &interior : *left_siblings) {
          if (interior != nullptr) {
            interior->latch_.unlock();
          }
        }
        for (const auto &interior : *right_siblings) {
          if (interior != nullptr) {
            interior->latch_.unlock();
          }
        }
        potential_changes->erase(potential_changes->begin(), potential_changes->end());
        left_siblings->erase(left_siblings->begin(), left_siblings->end());
        right_siblings->erase(right_siblings->begin(), right_siblings->end());
        indices->erase(indices->begin(), indices->end());
      }

      // Index in potential_changes.end() that will get you to next child
      i = FindKey(current, k) - 1;
      indices->push_back(i);
      left_siblings->push_back(left);
      right_siblings->push_back(right);
      potential_changes->push_back(current);

      // Only add the right sibling (potential right merge/borrow node) if leftmost node
      if (current->leaf_children_) {
        leaf = current->Leaf(i);
        leaf->latch_.lock();
        if (leaf->prev_ != nullptr) {
          leaf->prev_->latch_.lock();
        }
        if (leaf->next_ != nullptr) {
          leaf->next_->latch_.lock();
        }
      } else {
        if (i == 0 && left != nullptr) {
          left = left->Interior(left->filled_keys_ - 1);
          left->latch_.lock();
        } else if (i != 0) {
          left = current->Interior(i - 1);
          left->latch_.lock();
        }

        if (i == current->filled_keys_ - 1 && right != nullptr) {
          right = right->Interior(0);
          right->latch_.lock();
        } else if (i != current->filled_keys_ - 1) {
          right = current->Interior(i + 1);
          right->latch_.lock();
        }

        current = current->Interior(i);
        current->latch_.lock();
      }
    }

    return leaf;
  }

  /**
   * Removes a key/value pair from the overflow node.
   * @param overflow A pointer to the location of the first overflow node in the chain.
   * @param k The key to remove.
   * @param v A pointer to the value to remove, or nullptr if we should remove an arbitrary value.
   * @param result A pointer which gets filled with the value of the entry removed
   * @returns true if a value was removed, false otherwise.
   */
  bool RemoveFromOverflow(OverflowNode **overflow, KeyType k, ValueType *v, ValueType *result) {
    OverflowNode *current = *overflow;
    OverflowNode **prev_next = overflow;
    OverflowNode *prev = nullptr;
    bool deleted = false;
    uint32_t i = 0;

    while (current != nullptr) {
      for (i = 0; i < current->filled_keys_; i++) {
        if (KeyCmpEqual(k, current->keys_[i]) && (v == nullptr || value_eq_obj_(*v, current->values_[i]))) {
          *result = current->values_[i];
          deleted = true;
          break;
        }
      }
      if (i < current->filled_keys_) {
        break;
      }
      prev_next = &(current->next_);
      prev = current;
      current = current->next_;
    }

    if (deleted) {
      TERRIER_ASSERT(current != nullptr, "Should have a node");
      OverflowNode *hole = current;
      while (current->next_ != nullptr) {
        prev_next = &(current->next_);
        prev = current;
        current = current->next_;
      }

      // This shouldn't happen, but sometimes we don't delete nodes we should... this will clean that up
      if (current->filled_keys_ == 0) {
        *prev_next = current->next_;
        OverflowNode::Delete(current);
        TERRIER_ASSERT(prev != nullptr, "Must have at least one overflow");
        current = prev;
        // Don't have to update prev because this current will not get deleted
      }

      TERRIER_ASSERT(current->filled_keys_ > 0, "Should have keys in the node");
      hole->keys_[i] = current->keys_[current->filled_keys_ - 1];
      hole->values_[i] = current->values_[current->filled_keys_ - 1];
      current->filled_keys_--;
      if (current->filled_keys_ == 0) {
        OverflowNode::Delete(current);
        *prev_next = nullptr;
      }
    }

    return deleted;
  }

  /**
   * Merges two overflow nodes into one.
   * @param pointer The location of the primary OverflowNode chain to merge into.
   * @param other A pointer to the OverflowNode chain to merge in.
   */
  void MergeOverflow(OverflowNode **pointer, OverflowNode *other) {
    if (other == nullptr) return;  // No merging

    OverflowNode *last = nullptr;
    if (*pointer != nullptr) {
      last = *pointer;
      while (last->next_ != nullptr) {
        pointer = &last->next_;
        last = last->next_;
      }
    }

    // Join all full nodes
    *pointer = other;

    // Merge last two nodes
    while (other->next_ != nullptr) {
      other = other->next_;
    }

    if (last != nullptr) {
      while (last->filled_keys_ > 0 && other->filled_keys_ < OVERFLOW_SIZE) {
        other->keys_[other->filled_keys_] = last->keys_[last->filled_keys_ - 1];
        other->values_[other->filled_keys_] = last->values_[last->filled_keys_ - 1];
        other->filled_keys_++;
        last->filled_keys_--;
      }
      if (last->filled_keys_ == 0) {
        OverflowNode::Delete(last);
      } else {
        other->next_ = last;
      }
    }
  }

  /**
   * "Steals" elements with a specific key from one chain of OverflowNodes.
   * @param pointer The location of the primary OverflowNode chain to steal into.
   * @param other A pointer to the chain of OverflowNodes to steal from
   * @param key The key whose entries we should steal.
   */
  void StealOverflow(OverflowNode **pointer, OverflowNode *other, KeyType key) {
    if (other == nullptr) {
      return;  // Nothing to steal!
    }
    // Check if current node is full
    while (*pointer != nullptr && (*pointer)->filled_keys_ == OVERFLOW_SIZE) {
      pointer = &(*pointer)->next_;
    }
    OverflowNode *write_overflow = other;
    uint32_t write_index = 0;
    while (other != nullptr) {
      for (uint32_t i = 0; i < other->filled_keys_; i++) {
        if (KeyCmpEqual(key, other->keys_[i])) {
          // Move it to the new node
          if (*pointer == nullptr) {
            *pointer = OverflowNode::CreateNew();
          }
          TERRIER_ASSERT((*pointer)->filled_keys_ < OVERFLOW_SIZE, "Should be room in the overflow node");
          (*pointer)->keys_[(*pointer)->filled_keys_] = key;
          (*pointer)->values_[(*pointer)->filled_keys_] = other->values_[i];
          (*pointer)->filled_keys_++;
          if ((*pointer)->filled_keys_ == OVERFLOW_SIZE) {
            pointer = &(*pointer)->next_;
          }
        } else {
          // Keep it in the same node
          if (write_index == OVERFLOW_SIZE) {
            write_overflow = write_overflow->next_;
            write_index = 0;
          }

          write_overflow->keys_[write_index] = other->keys_[i];
          write_overflow->values_[write_index] = other->values_[i];
          write_index++;
        }
      }
      other = other->next_;
    }

    // Update filled_keys
    write_overflow->filled_keys_ = write_index;

    // Free non-used things
    OverflowNode *prev;
    OverflowNode *current = write_overflow->next_;
    write_overflow->next_ = nullptr;
    while (current != nullptr) {
      prev = current;
      current = current->next_;
      OverflowNode::Delete(prev);
    }
  }

  /**
   * Removes a key/value pair from a node.
   * @param node The node to remove from.
   * @param index The index to remove at.
   */
  template <typename Value>
  void RemoveFromNode(GenericNode<Value> *node, uint32_t index) {
    TERRIER_ASSERT(index < node->filled_keys_, "Not a valid index for removal!");
    for (uint32_t j = index + 1; j < node->filled_keys_; j++) {
      node->keys_[j - 1] = node->keys_[j];
      node->values_[j - 1] = node->values_[j];
    }
    node->filled_keys_--;
  }

  /**
   * Rebalances a given node, potentially taking keys/values from a node to the left or the right.
   * @param node The node which has not enough keys and needs more.
   * @param left The left sibling of the node, or nullptr for no left sibling.
   * @param right The right sibling of the node, or nullptr for no right sibling.
   * @param start The start index of keys in the node (1 or 0)
   * @param parent The parent of these nodes.
   * @param index The index of node in parent.
   * @param borrowed_key A pointer to be filled with the key if a borrow occurred.
   * @param borrowed_from A pointer to be filled with left or right if a borrow occurred from them.
   * @returns The primary node if "node" was deleted, or nullptr if a borrow occurred.
   */
  template <typename Value>
  GenericNode<Value> *Rebalance(GenericNode<Value> *node, GenericNode<Value> *left, GenericNode<Value> *right,
                                uint32_t start, InteriorNode *parent, uint32_t index, KeyType *borrowed_key,
                                GenericNode<Value> **borrowed_from) {
    TERRIER_ASSERT(node->filled_keys_ < MIN_CHILDREN, "Rebalance should be called on an unbalanced node!");

    if (index == 0) {  // 1 is the minimum index
      TERRIER_ASSERT(right != nullptr, "Can't have an empty right if we're deleting at index 0!");
      uint32_t size = right->filled_keys_;
      if (size > MIN_CHILDREN) {
        // Borrow case
        KeyType k_insert = start == 0 ? right->keys_[0] : parent->keys_[index + 1];
        parent->keys_[index + 1] = right->keys_[1];

        InsertIntoNode(node, node->filled_keys_, k_insert, right->values_[0]);
        RemoveFromNode(right, 0);

        if (borrowed_key != nullptr) {
          *borrowed_key = k_insert;
          *borrowed_from = right;
        }

        TERRIER_ASSERT(right->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
        TERRIER_ASSERT(node->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
        return nullptr;
      }
      // Merge case
      TERRIER_ASSERT(right->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
      // NB: using int64_t here in order to allow for natural semantics for a decreasing for loop
      for (int64_t i = right->filled_keys_ - 1; i >= 0; i--) {
        right->keys_[i + node->filled_keys_] = right->keys_[i];
        right->values_[i + node->filled_keys_] = right->values_[i];
      }

      if (start == 1) {
        right->keys_[node->filled_keys_] = parent->keys_[index + 1];
      }
      right->filled_keys_ += node->filled_keys_;
      TERRIER_ASSERT(right->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");

      for (uint32_t i = 0; i < node->filled_keys_; i++) {
        right->keys_[i] = node->keys_[i];
        right->values_[i] = node->values_[i];
      }

      // Cut the middleman out
      parent->values_[index].as_interior_ = reinterpret_cast<InteriorNode *>(right);
      TERRIER_ASSERT(right->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
      TERRIER_ASSERT(node->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
      return right;
    }
    TERRIER_ASSERT(left != nullptr, "Can't have an empty left if we're not deleting at index 0!");
    uint32_t size = left->filled_keys_;
    if (size > MIN_CHILDREN) {
      // Borrow case
      auto key = left->keys_[size - 1];
      if (start == 1) {
        // Borrow the goalpost from the parent
        node->keys_[0] = parent->keys_[index];
      }
      InsertIntoNode(node, 0, key, left->values_[size - 1]);
      RemoveFromNode(left, size - 1);
      parent->keys_[index] = key;

      if (borrowed_key != nullptr) {
        *borrowed_key = key;
        *borrowed_from = left;
      }

      TERRIER_ASSERT(left->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
      TERRIER_ASSERT(node->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
      return nullptr;
    }
    // Merge Case
    uint32_t orig_size = left->filled_keys_;
    for (uint32_t i = 0; i < node->filled_keys_; i++) {
      left->keys_[left->filled_keys_] = node->keys_[i];
      left->values_[left->filled_keys_] = node->values_[i];
      left->filled_keys_++;
    }

    if (start == 1) {
      left->keys_[orig_size] = parent->keys_[index];
    }

    TERRIER_ASSERT(left->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
    TERRIER_ASSERT(node->filled_keys_ <= NUM_CHILDREN, "Cant have too many filled Keys!");
    return left;
  }

  /**
   * Attempts to delete a key/value pair from a leaf, if it does not cause rebalancing.
   * @param leaf The leaf to delete from.
   * @param k The key to delete.
   * @param v The value to delete.
   * @param i The index of k/v in the leaf.
   * @returns nullopt if we need to rebalance, true if the deletion occurred with no rebalance,
   * and false if it did not occur.
   */
  std::optional<bool> TryDeleteLeaf(LeafNode *leaf, KeyType k, ValueType v, uint32_t i) {
    // Can't delete a key that doesn't exist in the tree
    if (i >= leaf->filled_keys_ || !KeyCmpEqual(k, leaf->keys_[i])) {
      return std::make_optional(false);
    }

    if (allow_duplicates_) {
      if (value_eq_obj_(v, leaf->values_[i])) {
        if (RemoveFromOverflow(&leaf->overflow_, k, nullptr, &leaf->values_[i])) {
          // We found a matching key/value pair!
          return std::make_optional(true);
        }
        // We did not find a matching key/value pair so we must fully delete this version
      } else {
        ValueType ignore;
        // Try to remove it from the overflow node.
        bool result = RemoveFromOverflow(&leaf->overflow_, k, &v, &ignore);
        return std::make_optional(result);
      }
    } else if (!value_eq_obj_(v, leaf->values_[i])) {
      return std::make_optional(false);  // Can't delete a key that doesn't match the value :(
    }

    if ((depth_ == 1 && root_->filled_keys_ == 1) || leaf->filled_keys_ > MIN_CHILDREN) {
      // No rebalancing needed
      RemoveFromNode(leaf, i);
      return std::make_optional(true);
    }

    return std::nullopt;
  }

  /**
   * Deletes the key/value pair from this tree.
   * @param k The key to delete.
   * @param v The value to delete.
   * @returns true if the delete was successful (i.e. k/v existed in the tree), false otherwise.
   */
  bool Delete(KeyType k, ValueType v) {
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif

    LeafNode *leaf = nullptr;
    uint32_t i;
    std::optional<bool> deleted;

    InteriorNode *save = root_;
    save->latch_.lock_shared();
    while (save != root_) {
      save->latch_.unlock_shared();
      save = root_;
      save->latch_.lock_shared();
    }

    leaf = TraverseTrack(root_, k, nullptr);

    i = FindKey(leaf, k);

    deleted = TryDeleteLeaf(leaf, k, v, i);
    leaf->latch_.unlock();
    if (deleted.has_value()) {
      // Successfully deleted without rebalancing!
#ifdef DEEP_DEBUG
      TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
      return deleted.value();
    }

    std::vector<InteriorNode *> potential_changes;
    std::vector<uint32_t> indices;
    std::vector<InteriorNode *> left_siblings;
    std::vector<InteriorNode *> right_siblings;

    save = root_;
    save->latch_.lock();
    while (save != root_) {
      save->latch_.unlock();
      save = root_;
      save->latch_.lock();
    }

    leaf = TraverseTrackWithSiblings(root_, k, &potential_changes, &indices, &left_siblings, &right_siblings);

    i = FindKey(leaf, k);

    auto unlock_all = [&]() {
      if (leaf->prev_ != nullptr) {
        leaf->prev_->latch_.unlock();
      }
      if (leaf->next_ != nullptr) {
        leaf->next_->latch_.unlock();
      }
      leaf->latch_.unlock();
      for (const auto &interior : potential_changes) {
        interior->latch_.unlock();
      }
      for (const auto &interior : left_siblings) {
        if (interior != nullptr) {
          interior->latch_.unlock();
        }
      }
      for (const auto &interior : right_siblings) {
        if (interior != nullptr) {
          interior->latch_.unlock();
        }
      }
    };

    deleted = TryDeleteLeaf(leaf, k, v, i);

    if (deleted.has_value()) {
      unlock_all();
#ifdef DEEP_DEBUG
      TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
      return deleted.value();
    }

    // At this point we are committed to deleting the key/value pair
    RemoveFromNode(leaf, i);

    LeafNode *borrowed_from;
    auto **parameter = reinterpret_cast<GenericNode<ValueType> **>(&borrowed_from);
    KeyType borrowed_key;
    auto *preserved = reinterpret_cast<LeafNode *>(Rebalance(
        leaf, leaf->prev_, leaf->next_, 0, potential_changes.back(), indices.back(), &borrowed_key, parameter));
    if (preserved == nullptr) {
      StealOverflow(&leaf->overflow_, borrowed_from->overflow_, borrowed_key);
      unlock_all();
#ifdef DEEP_DEBUG
      TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
      return true;  // No further rebalancing needed
    }

    TERRIER_ASSERT(preserved == leaf->prev_ || preserved == leaf->next_, "We always merge into leaf");

    // Merge our overflow into the preserved leaf
    MergeOverflow(&preserved->overflow_, leaf->overflow_);

    if (leaf->prev_ != nullptr) {
      leaf->prev_->next_ = leaf->next_;
      leaf->prev_->latch_.unlock();
    }
    if (leaf->next_ != nullptr) {
      leaf->next_->prev_ = leaf->prev_;
      leaf->next_->latch_.unlock();
    }

    if (preserved == leaf->prev_) {
      i = indices.back();
      indices.pop_back();
    } else {
      TERRIER_ASSERT(preserved == leaf->next_, "Only other option");
      i = indices.back() + 1;
      indices.pop_back();
      TERRIER_ASSERT(i < potential_changes.back()->filled_keys_, "i should always be a valid spot!");
    }
    LeafNode::Delete(leaf);

    InteriorNode *current;
    InteriorNode *left;
    InteriorNode *right;

    while (potential_changes.size() > 1) {
      current = potential_changes.back();
      left = left_siblings.back();
      right = right_siblings.back();
      potential_changes.pop_back();
      // We pop left/right at the bottom for easier coding

      TERRIER_ASSERT(left != nullptr || right != nullptr, "Should not have reached the root!");

      RemoveFromNode(current, i);
      GenericNode<Child> **unused = nullptr;
      auto *preserved_interior = reinterpret_cast<InteriorNode *>(
          Rebalance(current, left, right, 1, potential_changes.back(), indices.back(), nullptr, unused));
      if (preserved_interior == nullptr) {
        current->latch_.unlock();
        for (const auto &interior : potential_changes) {
          interior->latch_.unlock();
        }
        for (const auto &interior : left_siblings) {
          if (interior != nullptr) {
            interior->latch_.unlock();
          }
        }
        for (const auto &interior : right_siblings) {
          if (interior != nullptr) {
            interior->latch_.unlock();
          }
        }
#ifdef DEEP_DEBUG
        TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
        return true;  // No further rebalancing needed!
      }

      TERRIER_ASSERT(preserved_interior == left || preserved_interior == right, "We always merge into leaf");

      if (preserved_interior == left) {
        i = indices.back();
        indices.pop_back();
        TERRIER_ASSERT(potential_changes.back()->Interior(i) == current, "Child should still exist in parent");
      } else {
        TERRIER_ASSERT(preserved_interior == right, "Only other option");
        i = indices.back() + 1;
        indices.pop_back();
        TERRIER_ASSERT(i < potential_changes.back()->filled_keys_, "i should always be a valid spot!");
      }
      TERRIER_ASSERT(potential_changes.back()->latch_.holder_ == std::this_thread::get_id(),
                     "Lock should be held on parent for deletion");
      InteriorNode::Delete(current);

      // Unlock left/right latches
      if (left != nullptr) {
        left->latch_.unlock();
      }
      if (right != nullptr) {
        right->latch_.unlock();
      }

      left_siblings.pop_back();
      right_siblings.pop_back();
    }

    current = potential_changes.back();
    potential_changes.pop_back();
    TERRIER_ASSERT(potential_changes.empty(), "Should be done");

    left = left_siblings.back();
    right = right_siblings.back();
    if (left != nullptr) {
      left->latch_.unlock();
    }
    if (right != nullptr) {
      right->latch_.unlock();
    }
    RemoveFromNode(current, i);

    if (current->filled_keys_ >= MIN_CHILDREN) {
      current->latch_.unlock();
#ifdef DEEP_DEBUG
      TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
      return true;
    }

    TERRIER_ASSERT(current == root_, "Only root can be less than half full at this point");

    if (current->filled_keys_ == 1 && depth_ > 1) {
      --depth_;
      InteriorNode *interior = root_->Interior(0);
      for (uint32_t j = 0; j < interior->filled_keys_; j++) {
        root_->keys_[j] = interior->keys_[j];
        root_->values_[j] = interior->values_[j];
      }
      root_->filled_keys_ = interior->filled_keys_;
      root_->leaf_children_ = interior->leaf_children_;
      InteriorNode::Delete(interior);
      current->latch_.unlock();
    } else {
      current->latch_.unlock();
    }

#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "Deleting a key requires a valid B+tree");
#endif
    return true;
  }

  /**
   * Gets the heap usage of the entire tree, in bytes.
   */
  size_t GetHeapUsage() const {
#ifdef DEEP_DEBUG
    TERRIER_ASSERT(IsBplusTree(), "GetHeapUsage mulsled on a valid B+ Tree");
#endif
    return sizeof(BPlusTree) + root_->GetHeapUsage();
  }
};

template <typename KeyType, typename ValueType, typename KeyComparator, typename KeyEqualityChecker,
          typename KeyHashFunc, typename ValueEqualityChecker>
execution::sql::MemoryPool BPlusTree<KeyType, ValueType, KeyComparator, KeyEqualityChecker, KeyHashFunc,
                                     ValueEqualityChecker>::mem_pool{nullptr};

#undef CHECK
#undef CHECK_LT
#undef CHECK_LE
#undef CHECK_EQ
#undef DEEP_DEBUG
}  // namespace terrier::storage::index
