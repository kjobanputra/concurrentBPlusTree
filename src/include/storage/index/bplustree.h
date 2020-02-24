#pragma once

#include <functional>
#include <queue>

#include "common/spin_latch.h"

namespace terrier::storage::index {
// This macro is private to this file - see the #undef at the bottom
#define CHECK(x) do { if (!(x)) { return false; } } while(0)

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

    bool is_leaf_node(bool allow_duplicates) {
      CHECK(filled_keys_ >= 0 && filled_keys_ <= NUM_CHILDREN);

      for(int i = 1; i < filled_keys_; i++) {
          CHECK(keys_[i] < keys_[i-1]);
          CHECK(allow_duplicates || keys_[i] != keys_[i-1]);
      }
      return true;
    }
  };

  class InteriorNode {
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN - 1];
    bool leaf_children_;
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };

    bool is_interior_node(bool allow_duplicates,
                          InteriorNode *prev,
                          InteriorNode *next) {
      CHECK(filled_guide_posts_ > 0 && filled_guide_posts_ <= NUM_CHILDREN - 1);

      for(int i = 0; i <= filled_guide_posts_; i++) {
        if(leaf_children_) {
          CHECK(leaves_[i] != nullptr);
          CHECK(leaves_[i]->is_leaf_node(allow_duplicates));

          CHECK(i == 0 || (leaves_[i]->prev_ == leaves_[i-1] && leaves_[i-1]->next_ == leaves_[i]));

          auto leaf = leaves_[i];
          CHECK(i == 0 || (guide_posts_[i-1] <= leaf->keys_[0]));
          CHECK(i == filled_guide_posts_ || (leaf->keys_[leaf->filled_keys_-1] <= guide_posts_[i]));
        } else {
          CHECK(interiors_[i] != nullptr);
          CHECK(interiors_[i]->is_interior_node(true,
                                                   i == 0 ? prev : interiors_[i-1],
                                                   i == filled_guide_posts_ ? next : interiors_[i+1]));

          auto interior = interiors_[i];
          CHECK(i == 0 || (guide_posts_[i-1] <= interior->guide_posts_[0]));
          CHECK(i == filled_guide_posts_ || (interior->guide_posts_[interior->filled_guide_posts_-1] <= guide_posts_[i]));
        }
      }

      if(leaf_children_) {
        CHECK(prev == nullptr || prev->leaf_children_);
        CHECK((prev == nullptr && leaves_[0]->prev_ == nullptr) ||
              (prev != nullptr && leaves_[0]->prev_ == prev->leaves_[prev->filled_guide_posts_]));
        CHECK(prev == nullptr || (prev->guide_posts_[prev->filled_guide_posts_] <= leaves_[0]->keys_[0]));

        auto last_leaf = leaves_[filled_guide_posts_];
        CHECK(next == nullptr || next->leaf_children_);
        CHECK((next == nullptr && last_leaf->next_ == nullptr) ||
              (next != nullptr && last_leaf->next_ != next->leaves_[0]));

        CHECK(next == nullptr || (next->guide_posts_[0] <= last_leaf->keys_[last_leaf->filled_keys_-1]));
      } else {
        CHECK(prev == nullptr || (prev->guide_posts_[prev->filled_guide_posts_] <= interiors_[0]->guide_posts_[0]));
        auto last_interior = interiors_[filled_guide_posts_];
        CHECK(next == nullptr || (next->guide_posts_[0] <= last_interior->guide_posts_[last_interior->filled_guide_posts_]));
      }

      return true;
    }
  };

  InteriorNode *root_;
  uint32_t depth_;
  common::SpinLatch guard_;
  bool allow_duplicates_;

  bool is_bplus_tree() {
    if(root_->filled_guide_posts_ == 0) {
      // Root is not a valid interior node until we do a split!
      // However, this is only allowed to happen if the depth is truly 1
      return depth_ == 1 &&
             root_->leaf_children_ &&
             root_->leaves_[0]->is_leaf_node(allow_duplicates_) &&
             root_->leaves_[0]->next_ == nullptr &&
             root_->leaves[0]->prev_ == nullptr;
    }

    if(!root_->is_interior_node(allow_duplicates_, nullptr, nullptr)) {
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

 public:
  explicit BPlusTree(bool allow_duplicates): depth_(1), allow_duplicates_(allow_duplicates) {
    root_ = calloc(sizeof(class InteriorNode), 1);
    root_->leaf_children_ = true;
    root_->leaves[0] = calloc(sizeof(class LeafNode), 1);
  }

  bool insert(KeyType k, ValueType v) {
    TERRIER_ASSERT(is_bplus_tree(), "Insert must be called on a valid B+ Tree");
  }
};

#undef CHECK
}  // namespace terrier::storage::index
