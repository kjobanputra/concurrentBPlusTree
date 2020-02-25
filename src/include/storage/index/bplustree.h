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
    uint32_t filled_keys_;
    KeyType keys_[NUM_CHILDREN];
    ValueType values_[NUM_CHILDREN];
    LeafNode *prev_;
    LeafNode *next_;
    OverflowNode *overflow_;

    bool isLeafNode(bool allow_duplicates, bool is_root) {
      if(is_root) {
        CHECK(0 <= filled_keys_);
      } else {
        CHECK(MIN_CHILDREN <= filled_keys_);
      }
      CHECK(filled_keys_ <= NUM_CHILDREN);

      CHECK(allow_duplicates || overflow_ == nullptr);

      for(int i = 1; i < filled_keys_; i++) {
          CHECK(keys_[i] < keys_[i-1]);
          CHECK(keys_[i] != keys_[i-1]);
      }
      return true;
    }
  };

  class InteriorNode {
    uint32_t filled_guide_posts_;
    KeyType guide_posts_[NUM_CHILDREN - 1];
    union {
      InteriorNode *interiors_[NUM_CHILDREN];
      LeafNode *leaves_[NUM_CHILDREN];
    };
    bool leaf_children_;

    bool isInteriorNode(bool allow_duplicates,
                          InteriorNode *prev,
                          InteriorNode *next,
                          bool is_root) {
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
          CHECK(i == 0 || (guide_posts_[i-1] <= leaf->keys_[0]));
          CHECK(i == filled_guide_posts_ || (leaf->keys_[leaf->filled_keys_-1] <= guide_posts_[i]));
        } else {
          CHECK(interiors_[i] != nullptr);
          CHECK(interiors_[i]->isInteriorNode(true, i == 0 ? prev : interiors_[i - 1],
                                              i == filled_guide_posts_ ? next : interiors_[i + 1]));

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
  explicit BPlusTree(bool allow_duplicates): depth_(1), allow_duplicates_(allow_duplicates) {
    root_ = calloc(sizeof(class InteriorNode), 1);
    root_->leaf_children_ = true;
    root_->leaves[0] = calloc(sizeof(class LeafNode), 1);
  }

  bool Insert(KeyType k, ValueType v) {
    TERRIER_ASSERT(isBplusTree(), "Insert must be called on a valid B+ Tree");

    std::vector<InteriorNode *> potential_changes;
    potential_changes.reserve(depth_);
    InteriorNode *current = root_;
    LeafNode *leaf = nullptr;
    while(leaf == nullptr) {
      if(current->filled_guide_posts_ < NUM_CHILDREN) {
        potential_changes.erase(potential_changes.begin(), potential_changes.end());
      }
      potential_changes.push_back(current);

      uint32_t i = 0;
      while(i < current->filled_guide_posts_ && k >= current->guide_posts_[i]) {
        ++i;
      }

      if(current->leaf_children_) {
        leaf = current->leaves_[i];
        TERRIER_ASSERT(leaf != nullptr, "Leaves should not be null!!");
      } else {
        current = current->interiors_[i];
      }
    }

    if(leaf->filled_keys_ < NUM_CHILDREN) {
      potential_changes.erase(potential_changes.begin(), potential_changes.end());
    }

    uint32_t i = 0;
    while(i < leaf->filled_keys_ && k >= leaf->keys_[i]) {
      ++i;
    }

    if(!allow_duplicates_ && k == leaf->keys_[i]) {
      return false;
    } else if(k == leaf->keys_[i]) {
      OverflowNode *currentOverflow = leaf->overflow_;
      OverflowNode **prevOverflow = &(leaf->overflow_);
      while(currentOverflow != nullptr && currentOverflow->filled_keys_ >= OVERFLOW_SIZE) {
        prevOverflow = &(currentOverflow->next_);
        currentOverflow = currentOverflow->next_;
      }
      if(currentOverflow == nullptr) {
        currentOverflow = calloc(sizeof(OverflowNode), 1);
        *prevOverflow = currentOverflow;
      }

      currentOverflow->keys_[currentOverflow->filled_keys_] = k;
      currentOverflow->values_[currentOverflow->filled_keys_] = v;
      currentOverflow->filled_keys++;
      return true;
    }

    if(leaf->filled_keys_ == NUM_CHILDREN) {

      LeafNode *new_leaf = calloc(sizeof(LeafNode), 1);
      KeyType guide_post = split(reinterpret_cast<GenericNode<ValueType> *>(leaf),
                                 reinterpret_cast<GenericNode<ValueType> *>(new_leaf),
                                  i, k, v);

      TERRIER_ASSERT(leaf->isLeafNode(allow_duplicates_, false), "Old leaf was not preserved as leaf");
      TERRIER_ASSERT(new_leaf->isLeafNode(allow_duplicates_, false), "New Leaf not preserved as leaf");

      TERRIER_ASSERT(!potential_changes.empty(), "Potential changes should not be empty!!");
      InteriorNode *newInner = reinterpret_cast<LeafNode *>(new_leaf); // Pointers are the same size and we never dereference this!
      auto inner = --potential_changes.rend();
      for(; inner != potential_changes.rbegin(); --inner) {
        while(i < inner->filled_guide_posts_ && guide_post >= inner->guide_posts_[i]) {
          ++i;
        }

        InteriorNode *newestInner = calloc(sizeof(InteriorNode), 1);
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(newestInner),
                           i, guide_post, newInner);
        newInner = newestInner;
      }

      while(i < inner->filled_guide_posts_ && guide_post >= inner->guide_posts_[i]) {
        ++i;
      }
      if(inner->filled_guide_posts_ < NUM_CHILDREN - 1) {
        for(uint32_t j = inner->filled_keys_; j > i; --j) {
          inner->keys_[j] = inner->keys_[j-1];
          inner->interiors_[j] = inner->interiors_[j-1];
        }

        inner->keys_[i] = guide_post;
        inner->interiors_[i] = newInner;
      } else {
        TERRIER_ASSERT(inner == root_, "Top level inner potential change should not be full unless it is the root");

        InteriorNode *newestInner = calloc(sizeof(InteriorNode), 1);
        guide_post = split(reinterpret_cast<GenericNode<InteriorNode *>>(&(*inner)),
                           reinterpret_cast<GenericNode<InteriorNode *>>(newestInner),
                           i, guide_post, newInner);

        InteriorNode *newRoot = calloc(sizeof(InteriorNode), 1);
        newRoot->guide_posts_[0] = guide_post;
        newRoot->filled_guide_posts_ = 1;
        newRoot->interiors_[0] = inner;
        newRoot->interiors_[1] = newestInner;
        root_ = newRoot;
      }
    } else {
      for(uint32_t j = leaf->filled_keys_; j > i; --j) {
        leaf->keys_[j] = leaf->keys_[j-1];
        leaf->values_[j] = leaf->values_[j-1];
      }

      leaf->keys_[i] = k;
      leaf->values_[i] = v;
    }

    TERRIER_ASSERT(isBplusTree(), "End of insert should result in a valid B+ tree");
  }
};

#undef CHECK
}  // namespace terrier::storage::index
