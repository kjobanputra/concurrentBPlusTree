#pragma once

#include <functional>
#include <queue>
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
#define CHECK(x)    \
  do {              \
    if (!(x)) {     \
      return false; \
    }               \
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
 * Number of children an interior node is allowed to have, or
 * equivalently the number of keys a leaf node is allowed to store
 */
constexpr uint32_t NUM_CHILDREN = 5;

/**
 * Number of keys stored in an overflow node.
 */
constexpr uint32_t OVERFLOW_SIZE = 5;

/**
 * Minimum number of children an interior node is allowed to have
 * (equivalently, the minimum number of keys a leaf node is allowed to store)
 */
constexpr uint32_t MIN_CHILDREN = NUM_CHILDREN / 2 + NUM_CHILDREN % 2;

template <typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType>,
          typename KeyEqualityChecker = std::equal_to<KeyType>, typename KeyHashFunc = std::hash<KeyType>,
          typename ValueEqualityChecker = std::equal_to<ValueType>>
class BPlusTree {
 private:
  /**
   * This inner class defines what an overflow node looks like. An overflow node
   * is defined to only contain duplicate keys of the LeafNode it is attached to
   */
  class OverflowNode {
    uint32_t filled_keys_;
    KeyType keys_[OVERFLOW_SIZE];
    ValueType values_[OVERFLOW_SIZE];
    OverflowNode *next_;
  };

  /**
   * This inner class defines what a LeafNode is. The keys are stored in
   * sorted order with no duplicates within the keys_ array. If the
   * tree supports duplicate keys, then any duplicates will be stored
   * in the overflow node.
   */
  class LeafNode {
   private:
    uint32_t filled_keys_;
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;
    OverflowNode *overflow_;

    /**
     * Checks to make sure that the given node is a leaf.
     * @param allow_duplicates If duplicates are allowed in this leaf node
     * @param is_root If this leaf node is the only leaf node in the tree
     * @param parent The BPlusTree that this leaf node is a part of
     * @return whether this node is a valid leaf
     */
    bool IsLeafNode(bool allow_duplicates, bool is_root, BPlusTree *parent) {
      if (is_root) {
        CHECK_LE(0, filled_keys_);  // Root may even be empty
      } else {
        CHECK_LE(MIN_CHILDREN, filled_keys_);
      }
      CHECK_LE(filled_keys_, NUM_CHILDREN);

      CHECK(allow_duplicates || overflow_ == nullptr);

      for (int i = 1; i < filled_keys_; i++) {
        CHECK(parent->KeyCmpLess(keys_[i], keys_[i - 1]));
        CHECK(!parent->KeyCmpEqual(keys_[i], keys_[i - 1]));
      }

      // TODO(astanesc): Add checks to make sure overflow node is correct
      return true;
    }

    friend class BPlusTree;
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
  class InteriorNode {
   private:
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN];
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };
    bool leaf_children_;

    /**
     * Returns whether or not this node is an interior node
     *
     * @param allow_duplicates Does the BPlusTree allow duplicates?
     * @param prev the left sibling of this interior node, or {@code nullptr} if no left sibling exists
     * @param next the right sibling of this interior node, or {@code nullptr} if no right sibling exists
     */
    bool IsInteriorNode(bool allow_duplicates, InteriorNode *prev, InteriorNode *next, bool is_root,
                        BPlusTree *parent) {
      // Check bounds on number of filled guide posts
      if (is_root) {
        CHECK_LT(0, filled_guide_posts_);
      } else {
        CHECK_LT(MIN_CHILDREN, filled_guide_posts_);
      }
      CHECK_LE(filled_guide_posts_, NUM_CHILDREN - 1);

      // Check to make sure that every child is a valid node
      for (int i = 0; i <= filled_guide_posts_; i++) {
        if (leaf_children_) {
          CHECK(leaves_[i] != nullptr);
          CHECK(leaves_[i]->IsLeafNode(allow_duplicates));

          CHECK(i == 0 || (leaves_[i]->prev_ == leaves_[i - 1] && leaves_[i - 1]->next_ == leaves_[i]));
        } else {
          CHECK(interiors_[i] != nullptr);
          CHECK(interiors_[i]->IsInteriorNode(true, i == 0 ? prev : interiors_[i - 1],
                                              i == filled_guide_posts_ ? next : interiors_[i + 1]));
        }
      }

      // Make sure guide posts are in sorted order with no dupes
      for (int i = 1; i < filled_guide_posts_; i++) {
        CHECK_LT(guide_posts_[i - 1], guide_posts_[i]);
      }

      // Check to make sure that each child has keys that are in the correct range
      for (int i = 0; i <= filled_guide_posts_; i++) {
        if (leaf_children_) {
          auto leaf = leaves_[i];
          CHECK(i == 0 || parent->KeyCmpLessEqual(guide_posts_[i - 1], leaf->keys_[0]));
          CHECK(i == filled_guide_posts_ ||
                parent->KeyCmpLessEqual(leaf->keys_[leaf->filled_keys_ - 1], guide_posts_[i]));
        } else {
          auto interior = interiors_[i];
          CHECK(i == 0 || (parent->KeyCmpLessEqual(guide_posts_[i - 1], interior->guide_posts_[0])));
          CHECK(i == filled_guide_posts_ ||
                parent->KeyCmpLessEqual(interior->guide_posts_[interior->filled_guide_posts_ - 1], guide_posts_[i]));
        }
      }

      // Check to make sure that children on the edges (e.g. at indices 0 and filled_guide_posts_)
      // have keys in the correct ranges
      if (leaf_children_) {
        CHECK(prev == nullptr || prev->leaf_children_);
        CHECK((prev == nullptr && leaves_[0]->prev_ == nullptr) ||
              (prev != nullptr && leaves_[0]->prev_ == prev->leaves_[prev->filled_guide_posts_]));
        CHECK(prev == nullptr ||
              parent->KeyCmpLessEqual(prev->guide_posts_[prev->filled_guide_posts_], leaves_[0]->keys_[0]));

        auto last_leaf = leaves_[filled_guide_posts_];
        CHECK(next == nullptr || next->leaf_children_);
        CHECK((next == nullptr && last_leaf->next_ == nullptr) ||
              (next != nullptr && last_leaf->next_ != next->leaves_[0]));

        CHECK(next == nullptr ||
              parent->KeyCmpLessEqual(next->guide_posts_[0], last_leaf->keys_[last_leaf->filled_keys_ - 1]));
      } else {
        CHECK(prev == nullptr ||
              parent->KeyCmpLessEqual(prev->guide_posts_[prev->filled_guide_posts_], interiors_[0]->guide_posts_[0]));
        auto last_interior = interiors_[filled_guide_posts_];
        CHECK(next == nullptr ||
              parent->KeyCmpLessEqual(next->guide_posts_[0],
                                      last_interior->guide_posts_[last_interior->filled_guide_posts_]));
      }

      return true;
    }

    friend class BPlusTree;
  };

  InteriorNode *root_;
  uint32_t depth_;
  common::SpinLatch guard_;
  bool allow_duplicates_;
  KeyComparator key_cmp_obj_;
  KeyEqualityChecker key_eq_obj_;
  KeyHashFunc key_hash_obj_;
  ValueEqualityChecker value_eq_obj_;

  /**
   * Checks if this B+ Tree is truly a B+ Tree
   */
  bool IsBplusTree() {
    if (root_->filled_guide_posts_ == 0) {
      // Root is not a valid interior node until we do a split!
      // However, this is only allowed to happen if the depth is truly 1
      return depth_ == 1 && root_->leaf_children_ && root_->leaves_[0]->IsLeafNode(allow_duplicates_, true, this) &&
             root_->leaves_[0]->next_ == nullptr && root_->leaves[0]->prev_ == nullptr;
    }

    if (!root_->IsInteriorNode(allow_duplicates_, nullptr, nullptr, true, this)) {
      return false;
    }

    // Check to make sure that the depth of the tree is the same regardless of
    // which way you go down the tree
    std::queue<InteriorNode *> level;
    level.push(root_);
    for (int i = 1; i < depth_; i++) {
      std::queue<InteriorNode *> next_level;
      while (!level.empty()) {
        auto elem = level.pop();
        CHECK(!elem->leaf_children_);
        for (int j = 0; j <= elem->filled_guide_posts_; j++) {
          next_level.push(elem->interiors_[j]);
        }
      }
      level = std::move(next_level);
    }

    while (!level.empty()) {
      auto elem = level.pop();
      CHECK(elem->leaf_children);
    }
    return true;
  }

  /**
   * A generic node represents the beginning portion of both a LeafNode and an InteriorNode.
   * This allows us to write code that is agnostic to the specifics of a leaf node vs an interior node,
   * but wich only operates on keys and a generic value.
   */
  template <typename Value>
  struct GenericNode {
    uint32_t keycount_;
    KeyType keys_[NUM_CHILDREN];
    Value values_[NUM_CHILDREN];
  };

  /**
   * Make sure that doing a reinterpret-cast will actually lead to correct behavior!
   */
  static_assert(offsetof(LeafNode, filled_keys_) == offsetof(GenericNode<ValueType>, keycount_));
  static_assert(offsetof(LeafNode, keys_) == offsetof(GenericNode<ValueType>, keys_));
  static_assert(offsetof(LeafNode, values_) == offsetof(GenericNode<ValueType>, values_));

  static_assert(offsetof(InteriorNode, filled_guide_posts_) == offsetof(GenericNode<InteriorNode *>, keycount_));
  static_assert(offsetof(InteriorNode, guide_posts_) == offsetof(GenericNode<InteriorNode *>, keys_));
  static_assert(offsetof(InteriorNode, interiors_) == offsetof(GenericNode<InteriorNode *>, values_));

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
  const KeyType &SplitNode(GenericNode<Value> *from, GenericNode<Value> *to, uint32_t insertion_spot, const KeyType &k,
                           const Value &v) {
    TERRIER_ASSERT(from != nullptr && to != nullptr, "No nullptrs allowed here!");
    // NB: This assert cannot properly check for the node being completely full, since an interior node is
    //     full whenever the number of keys is NUM_CHILDREN - 1 and a leaf node when the number of keys is NUM_CHILDREN
    TERRIER_ASSERT(from->keycount_ >= NUM_CHILDREN - 1, "Split node should only be called on a full node!");
    TERRIER_ASSERT(to->keycount_ == 0, "Split node should be called to split into an empty node!");

    KeyType *guide_post = nullptr;
    if (insertion_spot <= NUM_CHILDREN / 2) {
      guide_post = &(from->keys_[NUM_CHILDREN / 2]);
    } else if (insertion_spot == NUM_CHILDREN / 2 + 1) {
      guide_post = &k;
    } else {
      guide_post = &(from->keys_[NUM_CHILDREN / 2 + 1]);
    }
    TERRIER_ASSERT(guide_post != nullptr, "Guide post should not be null!");

    if (insertion_spot <= NUM_CHILDREN / 2) {
      for (uint32_t j = NUM_CHILDREN / 2; j > insertion_spot; --j) {
        from->keys_[j] = from->keys_[j - 1];
        from->values_[j] = from->values_[j - 1];
      }
      from->keys_[insertion_spot] = k;
      from->values_[insertion_spot] = v;
      from->filled_keys_ = NUM_CHILDREN / 2 + 1;
    } else {
      for (uint32_t j = NUM_CHILDREN / 2; j < insertion_spot; ++j) {
        to->keys_[to->filled_keys_] = from->keys_[j];
        to->values_[to->filled_keys_] = from->values_[j];
        ++to->keycount_;
      }

      to->keys_[to->filled_keys_] = k;
      to->values_[to->filled_values_] = v;
      ++to->keycount_;

      for (uint32_t j = insertion_spot; j < NUM_CHILDREN; ++j) {
        to->keys_[to->filled_keys_] = from->keys_[j];
        to->values_[to->filled_values_] = from->values_[j];
        ++to->keycount_;
      }
    }
    return *guide_post;
  }

 public:
  explicit BPlusTree(bool allow_duplicates, KeyComparator key_cmp_obj = KeyComparator{},
                     KeyEqualityChecker key_eq_obj = KeyEqualityChecker{}, KeyHashFunc key_hash_obj = KeyHashFunc{},
                     ValueEqualityChecker value_eq_obj = ValueEqualityChecker{})
      : depth_(1),
        allow_duplicates_(allow_duplicates),
        key_cmp_obj_(key_cmp_obj),
        key_eq_obj_(key_eq_obj),
        key_hash_obj_(key_hash_obj),
        value_eq_obj_(value_eq_obj) {
    root_ = reinterpret_cast<InteriorNode *>(calloc(sizeof(class InteriorNode), 1));
    root_->leaf_children_ = true;
    root_->leaves_[0] = reinterpret_cast<LeafNode *>(calloc(sizeof(class LeafNode), 1));
  }

  /*
   * KeyCmpLess() - Compare two keys for "less than" relation
   *
   * If key1 < key2 return true
   * If not return false
   */
  inline bool KeyCmpLess(const KeyType &key1, const KeyType &key2) const { return key_cmp_obj_(key1, key2); }

  /*
   * KeyCmpEqual() - Compare a pair of keys for equality
   *
   * This functions compares keys for equality relation
   */
  inline bool KeyCmpEqual(const KeyType &key1, const KeyType &key2) const { return key_eq_obj_(key1, key2); }

  /*
   * KeyCmpGreaterEqual() - Compare a pair of keys for >= relation
   *
   * It negates result of keyCmpLess()
   */
  inline bool KeyCmpGreaterEqual(const KeyType &key1, const KeyType &key2) const { return !KeyCmpLess(key1, key2); }

  /*
   * KeyCmpGreater() - Compare a pair of keys for > relation
   *
   * It flips input for keyCmpLess()
   */
  inline bool KeyCmpGreater(const KeyType &key1, const KeyType &key2) const { return KeyCmpLess(key2, key1); }

  /*
   * KeyCmpLessEqual() - Compare a pair of keys for <= relation
   */
  inline bool KeyCmpLessEqual(const KeyType &key1, const KeyType &key2) const { return !KeyCmpGreater(key1, key2); }

  class LeafNodeIterator {
   private:
    friend class BPlusTree;

    LeafNode *current_;
    uint32_t index_;

    LeafNodeIterator(LeafNode *current, uint32_t index) : current_(current), index_(index) {}

   public:
    inline ValueType &Value() const {
      TERRIER_ASSERT(current_ != nullptr, "Value called on invalid iterator");
      return current_->values_[index_];
    }

    inline const KeyType &Key() const {
      TERRIER_ASSERT(current_ != nullptr, "Key called on invalid iterator");
      return current_->keys_[index_];
    }

    inline LeafNodeIterator &operator++() {
      index_++;
      if (index_ >= current_->filled_keys_) {
        current_ = current_->next_;
        index_ = 0;
      }
      return *this;
    }

    inline LeafNodeIterator operator++(int) {
      LeafNodeIterator copy = *this;
      this->operator++();
      return copy;
    }

    inline LeafNodeIterator &operator--() {
      if (index_ == 0) {
        current_ = current_->prev_;
        index_ = current_->filled_keys_ - 1;
      } else {
        index_--;
      }
      return *this;
    }

    inline LeafNodeIterator operator--(int) {
      LeafNodeIterator copy = *this;
      this->operator--();
      return copy;
    }

    inline bool operator==(const LeafNodeIterator &other) const {
      return other.current_ == this->current_ && other.index_ == this->index_;
    }

    inline bool operator!=(const LeafNodeIterator &other) const { return !this->operator==(other); }
  };

  /**
   * Gets an iterator pointing to the first key/value pair in the tree
   * @return the iterator
   */
  LeafNodeIterator begin() const {  // NOLINT for STL name compability
    InteriorNode *current = root_;
    while (!current->leaf_children_) {
      current = current->interiors_[0];
    }
    return {current->leaves_[0], 0};
  }

  /**
   * Gets an iterator pointing to the first key/value pair in the tree, whose key is >= the key passed in.
   * If the passed in key is larger than all keys in the tree, returns the end iterator.
   * @param key the Lower bound (inclusive) for the
   * @return the iterator
   */
  LeafNodeIterator begin(const KeyType &key) const {  // NOLINT for STL name compability
    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    while (leaf == nullptr) {
      uint32_t child = 0;
      while (child < current->filled_guide_posts_ && KeyCmpGreaterEqual(key, current->guide_posts_[child])) {
        child++;
      }
      if (current->leaf_children_) {
        leaf = current->leaves_[child];
      } else {
        current = current->interiors_[child];
      }
    }

    TERRIER_ASSERT(leaf != nullptr, "Leaf should be reached");
    for (uint32_t index = 0; index < leaf->filled_keys_; index++) {
      if (KeyCmpLessEqual(key, leaf->keys_[index])) {
        return {leaf, index};
      }
    }

    // Key exists in the next node over
    LeafNodeIterator ret = {leaf, leaf->filled_keys_ - 1};
    return ++ret;
  }

  inline const LeafNodeIterator end() const {  // NOLINT for STL name compability
    return {nullptr, 0};
  }

  bool Insert(KeyType k, ValueType v) {
    TERRIER_ASSERT(IsBplusTree(), "Insert must be called on a valid B+ Tree");

    std::vector<InteriorNode *> potential_changes;  // Mark the interior nodes that may be split
    potential_changes.reserve(depth_);

    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    while (leaf == nullptr) {
      // If we know we do not need to split (e.g. there are still open guide posts,
      // then we don't need to keep track of anything above us
      if (current->filled_guide_posts_ < NUM_CHILDREN - 1) {
        potential_changes.erase(potential_changes.begin(), potential_changes.end());
      }
      // However, the level below us may still split, so we may still change. Keep track of that
      potential_changes.push_back(current);

      uint32_t i = 0;
      while (i < current->filled_guide_posts_ && KeyCmpGreaterEqual(k, current->guide_posts_[i])) {
        ++i;
      }

      // If we've reached the leaf level, break out!
      if (current->leaf_children_) {
        leaf = current->leaves_[i];
        TERRIER_ASSERT(leaf != nullptr, "Leaves should not be null!!");
      } else {
        current = current->interiors_[i];
      }
    }

    // If the leaf definitely has enough room, then we don't need to split anything!
    if (leaf->filled_keys_ < NUM_CHILDREN) {
      potential_changes.erase(potential_changes.begin(), potential_changes.end());
    }

    uint32_t i = 0;
    while (i < leaf->filled_keys_ && KeyCmpGreater(k, leaf->keys_[i])) {
      ++i;
    }

    // If duplicates are not allowed, do not insert duplicates!
    if (!allow_duplicates_ && KeyCmpEqual(k, leaf->keys_[i])) {
      return false;
    }

    // If duplicates _are_ allowed, then insert into the overflow node, allocating a new one as necessary
    if (KeyCmpEqual(k, leaf->keys_[i])) {
      OverflowNode *current_overflow = leaf->overflow_;
      OverflowNode **prev_overflow = &(leaf->overflow_);
      while (current_overflow != nullptr && current_overflow->filled_keys_ >= OVERFLOW_SIZE) {
        prev_overflow = &(current_overflow->next_);
        current_overflow = current_overflow->next_;
      }
      if (current_overflow == nullptr) {
        current_overflow = reinterpret_cast<OverflowNode *>(calloc(sizeof(OverflowNode), 1));
        *prev_overflow = current_overflow;
      }

      current_overflow->keys_[current_overflow->filled_keys_] = k;
      current_overflow->values_[current_overflow->filled_keys_] = v;
      current_overflow->filled_keys++;
      return true;
    }

    // We need to split!
    if (leaf->filled_keys_ == NUM_CHILDREN) {
      auto *new_leaf = reinterpret_cast<LeafNode *>(calloc(sizeof(LeafNode), 1));
      KeyType guide_post = split(reinterpret_cast<GenericNode<ValueType> *>(leaf),
                                 reinterpret_cast<GenericNode<ValueType> *>(new_leaf), i, k, v);

      // These two variables represent where in leaf's overflow nodes we are reading from
      OverflowNode *read_overflow = leaf->overflow_;
      int read_id = 0;

      // These two variables represent where in leaf's overflow nodes we are writing to
      OverflowNode *write_overflow = leaf->overflow_;
      int write_id = 0;

      // These two variables represent where in new_leaf's overflow nodes we are writing to
      OverflowNode *to_overflow = nullptr;
      int to_id = 0;

      // Unlike in leaf, we might have to allocate new overflow leaves. This is where we should write the
      // pointer in order to ensure that the link is created - we do it this way to support lazy allocation
      OverflowNode **to_overflow_alloc = &(new_leaf->overflow_);

      // Skip over any overflow elements that don't need to be moved at all
      while (read_overflow != nullptr && KeyCmpLess(read_overflow->keys_[read_id], guide_post)) {
        read_id++;
        write_id++;
        if (read_id == OVERFLOW_SIZE || read_id == read_overflow->filled_keys_) {
          read_overflow = read_overflow->next_;
          read_id = 0;
          write_overflow = write_overflow->next_;
          write_id = 0;
        }
      }

      // For every other overflow element, move it either from leaf to leaf (in another slot) or from leaf to
      // new leaf (in a new slot) depending on where it should live
      while (read_overflow != nullptr) {
        if (KeyCmpLess(read_overflow->keys_[read_id], guide_post)) {
          // Move from leaf to leaf
          write_overflow->keys_[write_id] = read_overflow->keys_[read_id];
          write_overflow->values_[write_id] = read_overflow->values_[read_id];

          // Advance the write index
          write_id++;
          if (write_id == OVERFLOW_SIZE) {
            write_overflow = write_overflow->next_;
            write_id = 0;
          }
        } else {
          // Lazy allocation still means we have to allocate eventually, and now we know we need allocate
          if (to_overflow == nullptr) {
            to_overflow = reinterpret_cast<OverflowNode *>(calloc(sizeof(OverflowNode *), 1));
            *to_overflow_alloc = to_overflow;
            to_overflow_alloc = &(to_overflow->next_);
          }

          // Move from leaf to new_leaf.
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
        if (read_id == OVERFLOW_SIZE || read_id == read_overflow->filled_keys_) {
          read_overflow = read_overflow->next_;
          read_id = 0;
        }
      }

      // If we wrote any values to the new_leaf, we now might need to
      // free some overflow nodes!
      if (write_overflow != nullptr) {
        write_overflow->filled_keys_ = write_id;
        OverflowNode *temp = write_overflow;
        write_overflow = write_overflow->next_;
        temp->next_ = nullptr;
        while (write_overflow != nullptr) {
          temp = write_overflow->next_;
          free(write_overflow);
          write_overflow = temp;
        }
      }

      // If we did not completely fill up the new leaf overflow nodes, mark where we filled up to
      if (to_overflow != nullptr) {
        to_overflow->filled_keys_ = to_id;
      }

      TERRIER_ASSERT(leaf->IsLeafNode(allow_duplicates_, false), "Old leaf was not preserved as leaf");
      TERRIER_ASSERT(new_leaf->IsLeafNode(allow_duplicates_, false), "New Leaf not preserved as leaf");

      TERRIER_ASSERT(!potential_changes.empty(), "Potential changes should not be empty!!");
      TERRIER_ASSERT(
          potential_changes.front().filled_guide_posts_ < NUM_CHILDREN - 1 || potential_changes.front() == root_,
          "Potential changes should contain a front which has room for a new child, or the front should be root_");

      // to_insert represents the current node to insert into the parent.
      auto *to_insert =
          reinterpret_cast<InteriorNode *>(new_leaf);  // Pointers are the same size and we never dereference this!

      // inner represents the parent (where we insert)
      auto inner = --potential_changes.rend();
      for (; inner != potential_changes.rbegin(); ++inner) {
        // Find out where to insert this node
        while (i < inner->filled_guide_posts_ && KeyCmpGreater(guide_post, inner->guide_posts_[i])) {
          ++i;
        }

        TERRIER_ASSERT(!KeyCmpEqual(guide_post, inner->guide_posts_[i]), "We should not have duplicated guide posts!");

        // Perform the insert
        auto *new_inner = reinterpret_cast<InteriorNode *>(calloc(sizeof(InteriorNode), 1));
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(new_inner), i, guide_post, to_insert);
        to_insert = new_inner;
      }

      // Find out where to insert the final block
      while (i < inner->filled_guide_posts_ && KeyCmpGreater(guide_post, inner->guide_posts_[i])) {
        ++i;
      }

      TERRIER_ASSERT(!KeyCmpEqual(guide_post, inner->guide_posts_[i]), "We should not have duplicated guide posts!");

      // Can we insert here?
      if (inner->filled_guide_posts_ < NUM_CHILDREN - 1) {
        // Yes! Just insert!
        for (uint32_t j = inner->filled_keys_; j > i; --j) {
          inner->keys_[j] = inner->keys_[j - 1];
          inner->interiors_[j] = inner->interiors_[j - 1];
        }

        inner->keys_[i] = guide_post;
        inner->interiors_[i] = to_insert;
      } else {
        // No - we've reached the root and must split
        TERRIER_ASSERT(inner == root_, "Top level inner potential change should not be full unless it is the root");

        // Split the root
        auto *new_inner = reinterpret_cast<InteriorNode *>(calloc(sizeof(InteriorNode), 1));
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(new_inner), i, guide_post, to_insert);

        auto *new_root = reinterpret_cast<InteriorNode *>(calloc(sizeof(InteriorNode), 1));
        new_root->guide_posts_[0] = guide_post;
        new_root->filled_guide_posts_ = 1;
        new_root->interiors_[0] = inner;
        new_root->interiors_[1] = new_inner;
        root_ = new_root;
      }
    } else {
      for (uint32_t j = leaf->filled_keys_; j > i; --j) {
        leaf->keys_[j] = leaf->keys_[j - 1];
        leaf->values_[j] = leaf->values_[j - 1];
      }

      leaf->keys_[i] = k;
      leaf->values_[i] = v;
    }

    TERRIER_ASSERT(IsBplusTree(), "End of insert should result in a valid B+ tree");
  }
};

#undef CHECK
#undef CHECK_LT
#undef CHECK_LE
}  // namespace terrier::storage::index
