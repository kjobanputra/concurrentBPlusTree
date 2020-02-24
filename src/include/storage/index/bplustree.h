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
    uint32_t filled_keys_;
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;
  };

  class InteriorNode {
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN - 1];
    bool leaf_children_;
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };
  };

  InteriorNode *root_;
  uint32_t depth_;
  common::SpinLatch guard_;

 public:
  class LeafNodeIterator {
   private:
    friend class BPlusTree;

    LeafNode *current_;
    uint32_t index_;

    LeafNodeIterator(LeafNode *current, uint32_t index): current_(current), index_(index) {}

   public:
    inline ValueType &Value() const {
      TERRIER_ASSERT(!IsEnd(), "Value called on invalid iterator");
      return current_->values_[index_];
    }

    inline const KeyType &Key() const {
      TERRIER_ASSERT(!IsEnd(), "Key called on invalid iterator");
      return current_->keys[index_];
    }

    inline void operator++() {
      index_++;
      if (index_ >= current_->filled_keys_) {
        current_ = current_->next_;
        index_ = 0;
      }
    }

    inline LeafNodeIterator operator++(int) {
      LeafNodeIterator copy = *this;
      this->operator++();
      return copy;
    }

    inline void operator--() {
      if (index_ == 0) {
        current_ = current_->prev_;
        index_ = current_->filled_keys_ - 1;
      } else {
        index_--;
      }
    }

    inline LeafNodeIterator operator--(int) {
      LeafNodeIterator copy = *this;
      this->operator--();
      return copy;
    }

    inline bool operator==(const LeafNodeIterator &other) const {
      return other.current_ == this->current_ && other.index_ == this->index_;
    }

    inline bool IsEnd() const {
      return current_ == nullptr;
    }
  };

  /**
   * Gets an iterator pointing to the first key/value pair in the tree
   * @return the iterator
   */
  LeafNodeIterator Begin() const {
    InteriorNode *current = root_;
    for (uint32_t current_depth = 0; current_depth < depth_ - 1; current_depth++) {
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
  LeafNodeIterator Begin(KeyType key) const {
    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    for (uint32_t current_depth = 0; current_depth < depth_; current_depth++) {
      uint32_t child;
      for (child = 0; child < current->filled_guide_posts_
        && !KeyComparator(key, current->guide_posts_[child]); child++);
      if (current_depth == depth_ - 1) {
        leaf = current->leaves_[child];
      } else {
        current = current->interiors_[child];
      }
    }
    TERRIER_ASSERT(leaf != nullptr, "Leaf should be reached");
    for (uint32_t index = 0; index < leaf->filled_keys_; index++) {
      if (!KeyComparator(leaf->keys_[index], key)) {
        return {leaf, index};
      }
    }
    LeafNodeIterator ret = {leaf, leaf->filled_keys_ - 1};
    return ++ret;
  }
};

}  // namespace terrier::storage::index
