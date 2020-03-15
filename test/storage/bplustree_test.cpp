#include "test_util/test_harness.h"
#include "storage/index/bplustree.h"
#include "test_util/bwtree_test_util.h"
#include "test_util/multithread_test_util.h"
#include "test_util/test_harness.h"

constexpr uint32_t UNUSED_ATTRIBUTE NUM_LEAVES = 25;
constexpr uint32_t UNUSED_ATTRIBUTE NUM_KEYS_PER_LEAF = 4;
constexpr uint32_t NUM_INSERTIONS = 10000; // NUM_LEAVES*NUM_KEYS_PER_LEAF;

namespace terrier::storage::index {

struct BPlusTreeTests : public TerrierTest {
 public:
  const uint32_t num_threads_ = MultiThreadTestUtil::HardwareConcurrency() + (MultiThreadTestUtil::HardwareConcurrency() % 2);

  storage::index::BPlusTree<uint32_t, uint32_t> bplustree_;
  std::vector<uint32_t> key_vec_;
  BPlusTreeTests() : bplustree_(false) {}
  bool CheckDelete(uint32_t key, uint32_t val) {
    if (!bplustree_.IsBplusTree()) {
      return false;
    }

    auto iter = bplustree_.BeginLessEqual(key);
    if (iter == bplustree_.end() || iter.Key() != key) {
      iter.ReleaseLock();
      return true;
    }
    if (iter.Value() == val) {
      iter.ReleaseLock();
      return false;
    }
    for (const auto &value : iter) {
      if (value == val) {
        iter.ReleaseLock();
        return false;
      }
    }
    iter.ReleaseLock();
    return true;
  }

  uint32_t CountValues(uint32_t key) {
    auto iter = bplustree_.BeginLessEqual(key);
    if (iter == bplustree_.end() || iter.Key() != key) {
      iter.ReleaseLock();
      return 0;
    }
    uint32_t result = 1;
    for (const auto UNUSED_ATTRIBUTE &value : iter) {
      ++result;
    }
    iter.ReleaseLock();
    return result;
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
  return; // sppeeed
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

// NOLINTNEXTLINE
/**
 * Adapted from https://github.com/wangziqi2013/BwTree/blob/master/test/basic_test.cpp and
 * https://github.com/wangziqi2013/BwTree/blob/master/test/main.cpp
 *
 * Test Basic Insert/Delete/GetValue with different patterns and multi thread
 */
// NOLINTNEXTLINE
TEST_F(BPlusTreeTests, Interleaved) {
const uint32_t basic_test_key_num = 1024;

common::WorkerPool thread_pool(num_threads_, {});
thread_pool.Startup();

/*
 * InsertTest1() - Each threads inserts in its own consecutive key subspace
 *
 * The intervals of each thread does not intersect, therefore contention
 * is very small and this test is supposed to be very fast
 *
 * |---- thread 0 ----|---- thread 1----|----thread 2----| .... |---- thread n----|
 */
auto insert_test1 = [&](uint32_t id) {
  for (uint32_t i = id * basic_test_key_num; i < static_cast<uint32_t>(id + 1) * basic_test_key_num; i++) {
    bplustree_.Insert(i, i + 1);
    bplustree_.Insert(i, i + 2);
    bplustree_.Insert(i, i + 3);
    bplustree_.Insert(i, i + 4);
  }
};

/*
 * DeleteTest1() - Same pattern as InsertTest1()
 */
auto delete_test1 = [&](uint32_t id) {
  for (uint32_t i = id * basic_test_key_num; i < static_cast<uint32_t>(id + 1) * basic_test_key_num; i++) {
    bplustree_.Delete(i, i + 1);
    bplustree_.Delete(i, i + 2);
    bplustree_.Delete(i, i + 3);
    bplustree_.Delete(i, i + 4);
  }
};

/*
 * InsertTest2() - All threads collectively insert on the key space
 *
 * | t0 t1 t2 t3 .. tn | t0 t1 t2 t3 .. tn | t0 t1 .. | .. |  ... tn |
 *
 * This test is supposed to be slower since the contention is very high
 * between different threads
 */
auto insert_test2 = [&](uint32_t id) {
  for (uint32_t i = 0; i < basic_test_key_num; i++) {
    uint32_t key = num_threads_ * i + id;

    bplustree_.Insert(key, key + 1);
    bplustree_.Insert(key, key + 2);
    bplustree_.Insert(key, key + 3);
    bplustree_.Insert(key, key + 4);
  }
};

/*
 * DeleteTest2() - The same pattern as InsertTest2()
 */
auto delete_test2 = [&](uint32_t id) {
  for (uint32_t i = 0; i < basic_test_key_num; i++) {
    uint32_t key = num_threads_ * i + id;

    bplustree_.Delete(key, key + 1);
    bplustree_.Delete(key, key + 2);
    bplustree_.Delete(key, key + 3);
    bplustree_.Delete(key, key + 4);
  }
};

/*
 * DeleteGetValueTest() - Verifies all values have been deleted
 *
 * This function verifies on key_num * thread_num key space
 */
auto delete_get_value_test = [&]() {
  for (uint32_t i = 0; i < basic_test_key_num * num_threads_; i++) {
    EXPECT_EQ(CountValues(i), 0);
  }
};

/*
 * InsertGetValueTest() - Verifies all values have been inserted
 */
auto insert_get_value_test = [&]() {
  for (uint32_t i = 0; i < basic_test_key_num * num_threads_; i++) {
    EXPECT_EQ(CountValues(i), 4);
  }
};

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, insert_test2);

insert_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, delete_test1);

delete_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, insert_test1);

insert_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, delete_test2);

delete_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, insert_test1);

insert_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, delete_test1);

delete_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, insert_test2);

insert_get_value_test();

MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, delete_test2);

delete_get_value_test();
}

TEST_F(BPlusTreeTests, ScanDelete) {
  const uint32_t basic_test_key_num = 5 * 1024;
  const uint32_t num_insertions = 5;

  common::WorkerPool thread_pool(num_threads_, {});
  thread_pool.Startup();

  /*
   * insert_scan_delete() - Each threads inserts, scans, and deletes a bunch of keys
   */
  auto insert_scan_delete = [&](uint32_t id) {
    for (uint32_t run = 0; run < num_insertions; run++) {
      for (uint32_t i = 0; i < basic_test_key_num; i++) {
        bool result = bplustree_.Insert(i, id * basic_test_key_num + i, true);
        EXPECT_TRUE(result);
      }
      scan_start:
      auto scan_itr = bplustree_.begin();
      while (scan_itr != bplustree_.end()) {
        uint32_t key = scan_itr.Key();
        uint32_t target_value = id * basic_test_key_num + key;
        uint32_t count = 0;
        if (scan_itr.Value() == target_value) {
          count++;
        }
        for (auto &value : scan_itr) {
          if (value == target_value) {
            count++;
          }
        }
        EXPECT_EQ(count, 1);
        ++scan_itr;
        if(!scan_itr.Valid()) {
          goto scan_start;
        }
      }
      for (uint32_t i = 0; i < basic_test_key_num; i++) {
        bool result = bplustree_.Delete(i, id * basic_test_key_num + i);
        EXPECT_TRUE(result);
      }
      scan_start_2:
      scan_itr = bplustree_.begin();
      while (scan_itr != bplustree_.end()) {
        uint32_t key = scan_itr.Key();
        uint32_t target_value = id * basic_test_key_num + key;
        uint32_t count = 0;
        if (scan_itr.Value() == target_value) {
          count++;
        }
        for (auto &value : scan_itr) {
          if (value == target_value) {
            count++;
          }
        }
        EXPECT_EQ(count, 0);
        ++scan_itr;
        if(!scan_itr.Valid()) {
          goto scan_start_2;
        }
      }
    }
  };

  MultiThreadTestUtil::RunThreadsUntilFinish(&thread_pool, num_threads_, insert_scan_delete);
}

}  // namespace terrier::storage::index
