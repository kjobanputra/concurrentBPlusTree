#include "test_util/test_harness.h"
#include "storage/index/bplustree.h"

constexpr uint32_t UNUSED_ATTRIBUTE NUM_LEAVES = 25;
constexpr uint32_t UNUSED_ATTRIBUTE NUM_KEYS_PER_LEAF = 4;
constexpr uint32_t NUM_INSERTIONS = 100; // NUM_LEAVES*NUM_KEYS_PER_LEAF;

namespace terrier::storage::index {

struct BPlusTreeTests : public TerrierTest {
 public:
  storage::index::BPlusTree<uint32_t, uint32_t> bplustree_;
  std::vector<uint32_t> key_vec_;
  BPlusTreeTests() : bplustree_(false) {}
  bool CheckDelete(uint32_t key) {
    if (!bplustree_.IsBplusTree()) {
      return false;
    }

    auto iter = bplustree_.BeginLessEqual(key);
    return iter == bplustree_.end() || iter.Key() != key;
  }

};

// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, DeleteTest) {
  key_vec_.reserve(NUM_INSERTIONS);
  for (uint32_t i = 0; i < NUM_INSERTIONS; i++) {
    bplustree_.Insert(i, 1, false, [](uint32_t val){ return false; });
    key_vec_.push_back(i);
  }

  while (!key_vec_.empty()) {
    uint32_t del_index = rand() % key_vec_.size();
    auto key = key_vec_[del_index];
    bool deleted = bplustree_.Delete(key, 1);
    EXPECT_TRUE(deleted);
    key_vec_.erase(key_vec_.begin()+del_index);
    EXPECT_TRUE(CheckDelete(key));
  }
}
}  // namespace terrier::storage::index
