#ifndef DEFS_H
#define DEFS_H


#include <limits.h>
#include <cassert>
#include <cstdlib> // for rand()

#include <list>
#include <unordered_set>
#include <unordered_map>
#include <string>
#include <iostream>
#include <vector>

#ifndef FLT_MAX
#define FLT_MAX 10000000
#endif

#ifndef UINT_MAX
#define UINT_MAX ((uint)-1)
#endif


#define FAIL(x) {cout << x << endl; exit(1);}

// the LOWER the debuglevel, the MORE debug output!
#define DEBUGLEVEL 10

#if DEBUGLEVEL<1
#define DEBUG0(x) x;
#else
#define DEBUG0(x) ;
#endif

#if DEBUGLEVEL<2
#define DEBUG1(x) x;
#else
#define DEBUG1(x) ;
#endif

#if DEBUGLEVEL<3
#define DEBUG2(x) x;
#else
#define DEBUG2(x) ;
#endif

#if DEBUGLEVEL<4
#define DEBUG3(x) x;
#else
#define DEBUG3(x) ;
#endif

#if DEBUGLEVEL<5
#define DEBUG4(x) x;
#else
#define DEBUG4(x) ;
#endif

#if DEBUGLEVEL<6
#define DEBUG5(x) x;
#else
#define DEBUG5(x) ;
#endif


typedef unsigned int uint;


// generec merge of lists
template<class T>
std::list<T>& operator+=(std::list<T>& l1, const std::list<T>& l2){
  l1.insert(l1.end(), l2.begin(), l2.end());
  return l1;
}

// generic output for unordered sets of things, seperated by " "
template<class T, class Q>
std::ostream& operator<<(std::ostream& os, const std::unordered_set<T, Q>& l){
  os << '(';
  if(l.empty()) return os << ')';
  for(typename std::unordered_set<T, Q>::const_iterator i = l.begin(); i != l.end(); ++i)
    os << *i << ' ';
  return os << '\b' << ')';
}

// generic output for unordered maps of things, seperated by " "
template<class T, class Q, class R>
std::ostream& operator<<(std::ostream& os, const std::unordered_map<T, Q, R>& l){
  os << '(';
  if(l.empty()) return os << ')';
  for(typename std::unordered_map<T, Q, R>::const_iterator i = l.begin(); i != l.end(); ++i)
    os << *i << ' ';
  return os << '\b' << ')';
}


// generic output for unordered multisets of things, seperated by " "
template<class T>
std::ostream& operator<<(std::ostream& os, const std::unordered_multiset<T>& l){
  os << '(';
  if(l.empty()) return os << ')';
  for(typename std::unordered_multiset<T>::const_iterator i = l.begin(); i != l.end(); ++i)
    os << *i << ' ';
  return os << '\b' << ')';
}

// generic output for lists of things, seperated by " "
template<class T>
std::ostream& operator<<(std::ostream& os, const std::list<T>& l){
  os << '(';
  if(l.empty()) return os << ')';
  for(typename std::list<T>::const_iterator i = l.begin(); i != l.end(); ++i)
    os << *i << ' ';
  return os << '\b' << ')';
}
// generic output for pairs of things, seperated by " "
template<class T, class Q>
std::ostream& operator<<(std::ostream& os, const std::pair<T, Q>& l){
  return os << '(' << l.first << ',' << l.second<<')';
}

// O(n) erase by value for lists
template<class T>
inline void erase_all(std::list<T>& l, const T& x){
  for(typename std::list<T>::iterator i = l.begin(); i != l.end();) if(*i == x) i = l.erase(i); else ++i;
}

// O(n) find by value for lists
template<class T>
inline typename std::list<T>::const_iterator find_first(const std::list<T>& l, const T& x){
  for(typename std::list<T>::const_iterator i = l.begin(); i != l.end();) if(*i == x) return i; else ++i;
  return l.end();
}

// multiplication for pairs of integers to pairs of uint
template<class T>
inline std::pair<uint, uint> operator*(const std::pair<T, T>& P, const uint x){
  return std::make_pair(x * P.first, x * P.second);
}
// erase just the first occurance of something from a multiset
template<class T>
bool erase_once(std::unordered_multiset<T>& m, const T& x){
  typename std::unordered_multiset<T>::iterator i(m.find(x));
  if(i != m.end()) m.erase(i); else return false;
  return true;
}
// union of unordered sets
template<class T, class Q>
void unordered_set_union(std::unordered_set<T, Q>& s1, const std::unordered_set<T, Q>& s2){
  for(const T& x : s2)
    s1.insert(x);
}

// roll left/right, thanks to Didac Perez at StackOverflow
template<class T>
inline T ror(T x, const uint moves){ return (x >> moves) | (x << (sizeof(T)*8 - moves)); }
template<class T>
inline T rol(T x, const uint moves){ return (x << moves) | (x >> (sizeof(T)*8 - moves)); }

// first, abbreviate the merge of std::lists
template<class T>
std::list<T>& operator+=(std::list<T>& l, const std::string& s){
  l.push_back(s);
  return l;
}
template<class T>
std::list<T>& operator+=(std::list<T>& l, const uint& s){
  l.push_back(s);
  return l;
}
// generic output for std::lists of things, seperated by " "
template<class T>
std::ostream& operator<<(std::ostream& os, const std::vector<T>& l){
  os << "(";
  if(l.empty()) return os << ")";
  for(typename std::vector<T>::const_iterator i = l.begin(); i != l.end(); ++i)
    os << *i << " ";
  return os << "\b)";
}







// thanks Mr.Ree at StackOverflow
inline void split( const  std::string& theString,
            std::list<std::string>& theStringList,
            const  std::string& theDelimiter)
{
  assert(theDelimiter.size()>0);
  size_t  start = 0, end = 0;

  while(end != std::string::npos){
    end = theString.find( theDelimiter, start);
    // If at end, use length=maxLength.  Else use length=end-start.
    theStringList.push_back( theString.substr( start, (end == std::string::npos) ? std::string::npos : end - start));
    // If at end, use start=maxSize.  Else use start=end+delimiter.
    start = ( (end > (std::string::npos - theDelimiter.size())) ? std::string::npos : end + theDelimiter.size());
  }
}
// split without allowing double delimeters creating empty strings in the vector
inline void split_no_empty( const  std::string& theString,
                     std::list<std::string>& theStringList,
                     const  std::string& theDelimiter = " ")
{
  assert(theDelimiter.size()>0);
  size_t  start = 0, end = 0;

  while(end != std::string::npos){
    end = theString.find( theDelimiter, start);
    if(start != end){
      // If at end, use length=maxLength.  Else use length=end-start.
      theStringList.push_back( theString.substr( start, (end == std::string::npos) ? std::string::npos : end - start));
      // If at end, use start=maxSize.  Else use start=end+delimiter.
      start = ( (end > (std::string::npos - theDelimiter.size())) ? std::string::npos : end + theDelimiter.size());
    } else ++start;
  }
}

// randomized compare to use in conjunction with sort to get a random permutation
template<class T>
  struct rnd_compare{
    bool operator()(const T& a, const T& b){
    return (rand() % 1);
  }
};


#endif

