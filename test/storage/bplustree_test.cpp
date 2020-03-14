#include "test_util/test_harness.h"
#include "storage/index/bplustree.h"

constexpr uint32_t UNUSED_ATTRIBUTE NUM_LEAVES = 25;
constexpr uint32_t UNUSED_ATTRIBUTE NUM_KEYS_PER_LEAF = 4;
constexpr uint32_t NUM_INSERTIONS = 10000; // NUM_LEAVES*NUM_KEYS_PER_LEAF;

namespace terrier::storage::index {

struct BPlusTreeTests : public TerrierTest {
 public:
  storage::index::BPlusTree<uint32_t, uint32_t> bplustree_;
  std::vector<uint32_t> key_vec_;
  BPlusTreeTests() : bplustree_(false) {}
  bool CheckDelete(uint32_t key, uint32_t val) {
    if (!bplustree_.IsBplusTree()) {
      return false;
    }

    auto iter = bplustree_.BeginLessEqual(key);
    while(iter != bplustree_.end() && iter.Key() == key) {
      if(iter.Value() == val) {
        iter.ReleaseLock();
        return false;
      }
      ++iter;
    }
    iter.ReleaseLock();
    return true;
  }

};

// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, SimpleDeleteTest) {
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
    EXPECT_TRUE(CheckDelete(key, 1));
  }
}

// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, MixedDeleteInsertTest) {
  std::set<uint32_t> inserted_keys;

  bool insert = true;
  for(uint32_t i = 0; i < NUM_INSERTIONS * 5; i++) {
    if(rand() % 5 == 0) {
      // Switch phases
      insert = !insert;
    }

    if(inserted_keys.size() == 0) {
      // have to insert there's no key!
      insert = true;
    }

    if(inserted_keys.size() == NUM_INSERTIONS) {
      // Have to delete theres no space to insert!
      insert = false;
    }

    if(insert) {
      // Insert phase!
      uint32_t key;
      do {
        key = rand() % NUM_INSERTIONS;
      } while(inserted_keys.count(key) != 0);
      EXPECT_TRUE(bplustree_.Insert(key, 1, false, [](uint32_t val){ return false; }));
      inserted_keys.insert(key);
    } else {
      // Delete phase!
      uint32_t key = rand() % NUM_INSERTIONS;
      if(inserted_keys.count(key) == 1) {
        EXPECT_TRUE(bplustree_.Delete(key, 1));
        inserted_keys.erase(inserted_keys.find(key));
      } else {
        EXPECT_FALSE(bplustree_.Delete(key, 1));
      }
      EXPECT_TRUE(CheckDelete(key, 1));
    }
  }

}

// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, DuplicateTests) {
  std::map<uint32_t, std::map<uint32_t, uint32_t>> keyvals_;
  for(uint32_t i = 0; i < NUM_INSERTIONS; i++) {
    bplustree_.Insert(i, 1, true, [](uint32_t val){ return false; });
    keyvals_[i][1] = 1;
    while(rand() % 2 == 1) {
      uint32_t val = rand() % 3;
      bplustree_.Insert(i, val, true, [](uint32_t val) { return false; });
      keyvals_[i][val] ++;
    }
  }

  for(uint32_t i = 0; i < 2 * NUM_INSERTIONS; i++) {
    uint32_t k = rand() % NUM_INSERTIONS;
    uint32_t v = rand() % 3;

    bool deleted = bplustree_.Delete(k, v);
    if(keyvals_[k][v] > 0) {
      EXPECT_TRUE(deleted);
      keyvals_[k][v]--;
      if(keyvals_[k][v] == 0) {
        EXPECT_TRUE(CheckDelete(k, v));
      } else {
        EXPECT_TRUE(bplustree_.IsBplusTree());
      }
    } else {
      EXPECT_FALSE(deleted);
      EXPECT_TRUE(CheckDelete(k, v));
    }
  }
}
}  // namespace terrier::storage::index
