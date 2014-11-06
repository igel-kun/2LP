#ifndef REDUCTIONS_DEF_HPP
#define REDUCTIONS_DEF_HPP

namespace cr {
  enum reduction_type {TRR1, TRR2, TRR3, TRR4, TRR5, TRR6, PRR1, PRR2, PRR3, PRR4, PRR5, PRR6, PRR7, PRR8, Fin, YL };
}


inline ostream& operator<<(ostream& os, const cr::reduction_type& t){
  switch(t){
    case cr::TRR1: os <<"T1"; break;
    case cr::TRR2: os <<"T2"; break;
    case cr::TRR3: os <<"T3"; break;
    case cr::TRR4: os <<"T4"; break;
    case cr::TRR5: os <<"T5"; break;
    case cr::TRR6: os <<"T6"; break;
    case cr::PRR1: os <<"P1"; break;
    case cr::PRR2: os <<"P2"; break;
    case cr::PRR3: os <<"P3"; break;
    case cr::PRR4: os <<"P4"; break;
    case cr::PRR5: os <<"P5"; break;
    case cr::PRR6: os <<"P6"; break;
    case cr::PRR7: os <<"P7"; break;
    case cr::PRR8: os <<"P8"; break;
    case cr::Fin:  os <<"Fin"; break;
    case cr::YL:  os <<"YL"; break;
  }
  return os;
}



#endif
