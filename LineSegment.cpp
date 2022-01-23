#include"LineSegment.h"

    LineSegment :: LineSegment (v3 p0=v3(), v3 p1=v3(), int i=0) { 
      v[0] = p0; 
      v[1] = p1;
      index = i;
      vertical = false;  
      if ((v[1].x - v[0].x) != 0) {
        a = (v[1].y - v[0].y)/(v[1].x - v[0].x);
        b = (v[0].y - (a * v[0].x));
      }
      else {
        vertical = true;
      }
    }
    // bool LineSegment ::operator==(const LineSegment &ls) const { 
    //     return ((v[0] == ls.v[0]) && (v[1] == ls.v[1])); 
    // }

    // friend LineSegment :: ostream& operator<<(ostream& os, const LineSegment& ls) {
    //   os << "V0: (" << ls.v[0] << "); V1: (" << ls.v[1] << ")";
    //   return os;
    // }
