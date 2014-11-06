

namespace cr{
  enum cache_strategy_t {
    CACHE_LFU,   // least frequently
    CACHE_LRU,   // least recently
    CACHE_MRU,   // most recently
  };
  struct cache_opts {
    uint size;
    cache_strategy_t strategy;
  };

  typedef map<const graph,const solution_t> solution_cache_t;

  // the global solution cache
  solution_cache_t solution_cache;

  // query the cache, returns the cached solution for this graph or solution_t() if none was found
  solution_t query_cache(const graph& g);

  // insert a solution for a graph into the cache
  void insert_into_cache(const graph& g, const solution_t sol);
};
