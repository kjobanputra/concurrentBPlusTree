#pragma once

#include <functional>
#include <queue>
#include <stack>

#include "common/spin_latch.h"

namespace terrier::storage::index {
// This macro is private to this file - see the #undef at the bottom
#define CHECK(x) do { if (!(x)) { return false; } } while(0)

constexpr uint32_t NUM_CHILDREN = 5;
constexpr uint32_t OVERFLOW_SIZE = 5;
constexpr uint32_t MIN_CHILDREN = NUM_CHILDREN / 2 + NUM_CHILDREN % 2;

template <typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType>,
          typename KeyEqualityChecker = std::equal_to<KeyType>, typename KeyHashFunc = std::hash<KeyType>,
          typename ValueEqualityChecker = std::equal_to<ValueType>>
class BPlusTree {
 private:
  class OverflowNode {
    uint32_t filled_keys_;
    KeyType keys_[OVERFLOW_SIZE];
    ValueType values_[OVERFLOW_SIZE];
    OverflowNode *next_;
  };

  class LeafNode {
   private:
    uint32_t filled_keys_;
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;
    OverflowNode *overflow_;

    bool isLeafNode(bool allow_duplicates, bool is_root, BPlusTree *parent) {
      if(is_root) {
        CHECK(0 <= filled_keys_);
      } else {
        CHECK(MIN_CHILDREN <= filled_keys_);
      }
      CHECK(filled_keys_ <= NUM_CHILDREN);

      CHECK(allow_duplicates || overflow_ == nullptr);

      for(int i = 1; i < filled_keys_; i++) {
          CHECK(parent->KeyCmpLess(keys_[i], keys_[i-1]));
          CHECK(!parent->KeyCmpEqual(keys_[i], keys_[i-1]));
      }
      return true;
    }

    friend class BPlusTree;
  };

  // NB: We allocate one more guide post than necessary in order to allow
  //     InteriorNode and LeafNode to have the same start (e.g. GenericNode)
  class InteriorNode {
   private:
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN];
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };
    bool leaf_children_;

    bool isInteriorNode(bool allow_duplicates,
                          InteriorNode *prev,
                          InteriorNode *next,
                          bool is_root,
                          BPlusTree *parent) {
      if(is_root) {
        CHECK(0 < filled_guide_posts_);
      } else {
        CHECK(MIN_CHILDREN < filled_guide_posts_);
      }

      CHECK(filled_guide_posts_ <= NUM_CHILDREN - 1);

      for(int i = 0; i <= filled_guide_posts_; i++) {
        if(leaf_children_) {
          CHECK(leaves_[i] != nullptr);
          CHECK(leaves_[i]->isLeafNode(allow_duplicates));

          CHECK(i == 0 || (leaves_[i]->prev_ == leaves_[i-1] && leaves_[i-1]->next_ == leaves_[i]));

          auto leaf = leaves_[i];
          CHECK(i == 0 || parent->KeyCmpLessEqual(guide_posts_[i-1], leaf->keys_[0]));
          CHECK(i == filled_guide_posts_ || parent->KeyCmpLessEqual(leaf->keys_[leaf->filled_keys_-1], guide_posts_[i]));
        } else {
          CHECK(interiors_[i] != nullptr);
          CHECK(interiors_[i]->isInteriorNode(true, i == 0 ? prev : interiors_[i - 1],
                                              i == filled_guide_posts_ ? next : interiors_[i + 1]));

          auto interior = interiors_[i];
          CHECK(i == 0 || (parent->KeyCmpLessEqual(guide_posts_[i-1], interior->guide_posts_[0])));
          CHECK(i == filled_guide_posts_ || parent->KeyCmpLessEqual(interior->guide_posts_[interior->filled_guide_posts_-1], guide_posts_[i]));
        }
      }

      if(leaf_children_) {
        CHECK(prev == nullptr || prev->leaf_children_);
        CHECK((prev == nullptr && leaves_[0]->prev_ == nullptr) ||
              (prev != nullptr && leaves_[0]->prev_ == prev->leaves_[prev->filled_guide_posts_]));
        CHECK(prev == nullptr || parent->KeyCmpLessEqual(prev->guide_posts_[prev->filled_guide_posts_], leaves_[0]->keys_[0]));

        auto last_leaf = leaves_[filled_guide_posts_];
        CHECK(next == nullptr || next->leaf_children_);
        CHECK((next == nullptr && last_leaf->next_ == nullptr) ||
              (next != nullptr && last_leaf->next_ != next->leaves_[0]));

        CHECK(next == nullptr || parent->KeyCmpLessEqual(next->guide_posts_[0], last_leaf->keys_[last_leaf->filled_keys_-1]));
      } else {
        CHECK(prev == nullptr || parent->KeyCmpLessEqual(prev->guide_posts_[prev->filled_guide_posts_], interiors_[0]->guide_posts_[0]));
        auto last_interior = interiors_[filled_guide_posts_];
        CHECK(next == nullptr || parent->KeyCmpLessEqual(next->guide_posts_[0], last_interior->guide_posts_[last_interior->filled_guide_posts_]));
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


  bool isBplusTree() {
    if(root_->filled_guide_posts_ == 0) {
      // Root is not a valid interior node until we do a split!
      // However, this is only allowed to happen if the depth is truly 1
      return depth_ == 1 &&
             root_->leaf_children_ && root_->leaves_[0]->isLeafNode(allow_duplicates_, true) &&
             root_->leaves_[0]->next_ == nullptr &&
             root_->leaves[0]->prev_ == nullptr;
    }

    if(!root_->isInteriorNode(allow_duplicates_, nullptr, nullptr, true)) {
      return false;
    }

    std::queue<InteriorNode *> level;
    level.push(root_);
    for(int i = 1; i < depth_; i++) {
      std::queue<InteriorNode *> next_level;
      while(!level.empty()) {
        auto elem = level.pop();
        CHECK(!elem->leaf_children_);
        for(int j = 0; j <= elem->filled_guide_posts_; j++) {
          next_level.push(elem->interiors_[j]);
        }
      }
      level = std::move(next_level);
    }

    while(!level.empty()) {
      auto elem = level.pop();
      CHECK(elem->leaf_children);
    }
    return true;
  }

  template <typename Value>
  struct GenericNode {
    uint32_t keycount_;
    KeyType keys_[NUM_CHILDREN];
    Value values_[NUM_CHILDREN];
  };

  static_assert(offsetof(LeafNode, filled_keys_) == offsetof(GenericNode<ValueType>, keycount_));
  static_assert(offsetof(LeafNode, keys_) == offsetof(GenericNode<ValueType>, keys_));
  static_assert(offsetof(LeafNode, values_) == offsetof(GenericNode<ValueType>, values_));

  static_assert(offsetof(InteriorNode, filled_guide_posts_) == offsetof(GenericNode<InteriorNode *>, keycount_));
  static_assert(offsetof(InteriorNode, guide_posts_) == offsetof(GenericNode<InteriorNode *>, keys_));
  static_assert(offsetof(InteriorNode, interiors_) == offsetof(GenericNode<InteriorNode *>, values_));

  template <typename Value>
  const KeyType &splitNode(GenericNode<Value> *from,
                           GenericNode<Value> *to,
                           uint32_t insertion_spot,
                           const KeyType &k,
                           const Value &v) {
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
  explicit BPlusTree(bool allow_duplicates,
                     KeyComparator key_cmp_obj = KeyComparator{},
                     KeyEqualityChecker key_eq_obj = KeyEqualityChecker{}, KeyHashFunc key_hash_obj = KeyHashFunc{},
                     ValueEqualityChecker value_eq_obj = ValueEqualityChecker{})
      : allow_duplicates_(allow_duplicates),
        key_cmp_obj_(key_cmp_obj), key_eq_obj_(key_eq_obj),
        key_hash_obj_(key_hash_obj), value_eq_obj_(value_eq_obj) {
    root_ = calloc(sizeof(class InteriorNode), 1);
    root_->leaf_children_ = true;
    root_->leaves[0] = calloc(sizeof(class LeafNode), 1);
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

    LeafNodeIterator(LeafNode *current, uint32_t index): current_(current), index_(index) {}

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

    inline bool operator!=(const LeafNodeIterator &other) const {
      return !this->operator==(other);
    }
  };

  /**
   * Gets an iterator pointing to the first key/value pair in the tree
   * @return the iterator
   */
  LeafNodeIterator begin() const { // NOLINT for STL name compability
    InteriorNode *current = root_;
    while(!current->leaf_children_) {
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
  LeafNodeIterator begin(const KeyType &key) const { // NOLINT for STL name compability
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

  inline const LeafNodeIterator end() const { // NOLINT for STL name compability
    return {nullptr, 0};
  }

  bool Insert(KeyType k, ValueType v) {
    TERRIER_ASSERT(isBplusTree(), "Insert must be called on a valid B+ Tree");

    std::vector<InteriorNode *> potential_changes;
    potential_changes.reserve(depth_);
    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    while (leaf == nullptr) {
      if (current->filled_guide_posts_ < NUM_CHILDREN) {
        potential_changes.erase(potential_changes.begin(), potential_changes.end());
      }
      potential_changes.push_back(current);

      uint32_t i = 0;
      while (i < current->filled_guide_posts_ && KeyCmpGreaterEqual(k, current->guide_posts_[i])) {
        ++i;
      }

      if (current->leaf_children_) {
        leaf = current->leaves_[i];
        TERRIER_ASSERT(leaf != nullptr, "Leaves should not be null!!");
      } else {
        current = current->interiors_[i];
      }
    }

    if (leaf->filled_keys_ < NUM_CHILDREN) {
      potential_changes.erase(potential_changes.begin(), potential_changes.end());
    }

    uint32_t i = 0;
    while (i < leaf->filled_keys_ && KeyCmpGreaterEqual(k, leaf->keys_[i])) {
      ++i;
    }

    if (!allow_duplicates_ && KeyCmpEqual(k, leaf->keys_[i])) {
      return false;
    } else if (KeyCmpEqual(k, leaf->keys_[i])) {
      OverflowNode *currentOverflow = leaf->overflow_;
      OverflowNode **prevOverflow = &(leaf->overflow_);
      while (currentOverflow != nullptr && currentOverflow->filled_keys_ >= OVERFLOW_SIZE) {
        prevOverflow = &(currentOverflow->next_);
        currentOverflow = currentOverflow->next_;
      }
      if (currentOverflow == nullptr) {
        currentOverflow = calloc(sizeof(OverflowNode), 1);
        *prevOverflow = currentOverflow;
      }

      currentOverflow->keys_[currentOverflow->filled_keys_] = k;
      currentOverflow->values_[currentOverflow->filled_keys_] = v;
      currentOverflow->filled_keys++;
      return true;
    }

    if (leaf->filled_keys_ == NUM_CHILDREN) {
      LeafNode *new_leaf = calloc(sizeof(LeafNode), 1);
      KeyType guide_post = split(reinterpret_cast<GenericNode<ValueType> *>(leaf),
                                 reinterpret_cast<GenericNode<ValueType> *>(new_leaf), i, k, v);

      TERRIER_ASSERT(leaf->isLeafNode(allow_duplicates_, false), "Old leaf was not preserved as leaf");
      TERRIER_ASSERT(new_leaf->isLeafNode(allow_duplicates_, false), "New Leaf not preserved as leaf");

      TERRIER_ASSERT(!potential_changes.empty(), "Potential changes should not be empty!!");
      InteriorNode *newInner =
          reinterpret_cast<LeafNode *>(new_leaf);  // Pointers are the same size and we never dereference this!
      auto inner = --potential_changes.rend();
      for (; inner != potential_changes.rbegin(); --inner) {
        while (i < inner->filled_guide_posts_ && KeyCmpGreaterEqual(guide_post, inner->guide_posts_[i])) {
          ++i;
        }

        InteriorNode *newestInner = calloc(sizeof(InteriorNode), 1);
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(newestInner), i, guide_post, newInner);
        newInner = newestInner;
      }

      while (i < inner->filled_guide_posts_ && KeyCmpGreaterEqual(guide_post, inner->guide_posts_[i])) {
        ++i;
      }
      if (inner->filled_guide_posts_ < NUM_CHILDREN - 1) {
        for (uint32_t j = inner->filled_keys_; j > i; --j) {
          inner->keys_[j] = inner->keys_[j - 1];
          inner->interiors_[j] = inner->interiors_[j - 1];
        }

        inner->keys_[i] = guide_post;
        inner->interiors_[i] = newInner;
      } else {
        TERRIER_ASSERT(inner == root_, "Top level inner potential change should not be full unless it is the root");

        InteriorNode *newestInner = calloc(sizeof(InteriorNode), 1);
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(newestInner), i, guide_post, newInner);

        InteriorNode *newRoot = calloc(sizeof(InteriorNode), 1);
        newRoot->guide_posts_[0] = guide_post;
        newRoot->filled_guide_posts_ = 1;
        newRoot->interiors_[0] = inner;
        newRoot->interiors_[1] = newestInner;
        root_ = newRoot;
      }
    } else {
      for (uint32_t j = leaf->filled_keys_; j > i; --j) {
        leaf->keys_[j] = leaf->keys_[j - 1];
        leaf->values_[j] = leaf->values_[j - 1];
      }

      leaf->keys_[i] = k;
      leaf->values_[i] = v;
    }

    TERRIER_ASSERT(isBplusTree(), "End of insert should result in a valid B+ tree");
  }
};

#undef CHECK
}  // namespace terrier::storage::index
