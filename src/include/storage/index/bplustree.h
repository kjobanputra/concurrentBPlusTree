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
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;
  };
  class InteriorNode {
    KeyType guideposts_[NUM_CHILDREN - 1];
    bool leaf_children_;
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };
  };

  InteriorNode *root_;
  uint32_t depth_;
  common::SpinLatch guard_;
};

}  // namespace terrier::storage::index
