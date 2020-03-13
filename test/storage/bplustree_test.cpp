#include "test_util/test_harness.h"
#include "storage/index/bplustree.h"

constexpr uint32_t NUM_LEAVES = 25;
constexpr uint32_t NUM_KEYS_PER_LEAF = 4;
constexpr uint32_t NUM_INSERTIONS = NUM_LEAVES*NUM_KEYS_PER_LEAF;

namespace terrier::storage::index {

struct BPlusTreeTests : public TerrierTest {
 public:
  storage::index::BPlusTree<uint32_t, uint32_t> *bplustree_;
  std::vector<uint32_t> keyVec_;
  BPlusTreeTests() : bplustree_(false) {}
  bool CheckDelete(uint32_t del_index) {
    bool is_b_plus_tree = bplustree_->IsBplusTree();
    /*
    for (auto it = keyVec_.begin(); it != keyVec_.end(); it++) {
      auto current = bplustree_->root_;

    }
    */

    return is_b_plus_tree;
  }

};

// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, DeleteTest) {
  keyVec_.reserve(NUM_INSERTIONS);
  for (uint32_t i = 0; i < NUM_INSERTIONS; i++) {
    bplustree_->Insert(i, 1, false, [](uint32_t val){ return false; });
    keyVec_.push_back(i);
  }

  while (!keyVec_.empty()) {
    uint32_t del_index = rand() % keyVec_.size();
    bplustree_->Delete(keyVec_[del_index]);
    keyVec_.erase(keyVec_.begin()+del_index);
    EXPECT_TRUE(CheckDelete(del_index));
  }
}
}  // namespace terrier::storage::index
