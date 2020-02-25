#pragma once

#include <functional>
#include "common/spin_latch.h"

namespace terrier::storage::index {

constexpr uint32_t NUM_CHILDREN = 5;

template <typename KeyType, typename ValueType, typename KeyComparator = std::less<KeyType>,
          typename KeyEqualityChecker = std::equal_to<KeyType>, typename KeyHashFunc = std::hash<KeyType>,
          typename ValueEqualityChecker = std::equal_to<ValueType>>
class BPlusTree {
 private:
  class LeafNode {
   private:
    uint32_t filled_keys_;
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;

    friend class BPlusTree;
  };

  class InteriorNode {
   private:
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN - 1];
    bool leaf_children_;
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };

    friend class BPlusTree;
  };

  InteriorNode *root_;
  uint32_t depth_;
  common::SpinLatch guard_;

  KeyComparator key_cmp_obj_;
  KeyEqualityChecker key_eq_obj_;
  KeyHashFunc key_hash_obj_;
  ValueEqualityChecker value_eq_obj_;

 public:
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
    do {
      current = current->interiors_[0];
    } while(!current->leaf_children_);
    return {current->leaves_[0], 0};
  }

  /**
   * Gets an iterator pointing to the first key/value pair in the tree, whose key is >= the key passed in.
   * If the passed in key is larger than all keys in the tree, returns the end iterator.
    * @param key the Lower bound (inclusive) for the
   * @return the iterator
   */
  LeafNodeIterator begin(KeyType key) const { // NOLINT for STL name compability
    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    do {
      uint32_t child = 0;
      while (child < current->filled_guide_posts_ && KeyCmpGreaterEqual(key, current->guide_posts_[child])) {
        child++;
      }
      if (current->leaf_children_) {
        leaf = current->leaves_[child];
        break;
      }
      current = current->interiors_[child];
    } while(true);

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

  explicit BPlusTree(KeyComparator key_cmp_obj = KeyComparator{},
         KeyEqualityChecker key_eq_obj = KeyEqualityChecker{}, KeyHashFunc key_hash_obj = KeyHashFunc{},
         ValueEqualityChecker value_eq_obj = ValueEqualityChecker{})
      : key_cmp_obj_(key_cmp_obj), key_eq_obj_(key_eq_obj),
      key_hash_obj_(key_hash_obj), value_eq_obj_(value_eq_obj) {
    //TODO: More stuff
  }
};

}  // namespace terrier::storage::index
