function [RotMat] = Quat2RotMat(QVec)

  RotMat = zeros(3) ;
  
  a = QVec(1) ;
  b = QVec(2) ;
  c = QVec(3) ;
  d = QVec(4) ;

  RotMat(1,1) = a^2 + b^2 - c^2 - d^2 ;
  RotMat(1,2) = 2*b*c + 2*a*d ;
  RotMat(1,3) = 2*b*d - 2*c*a ;
  RotMat(2,1) = 2*b*c - 2*a*d ;
  RotMat(2,2) = a^2 - b^2 + c^2 - d^2 ;
  RotMat(2,3) = 2*c*d + 2*a*b ;
  RotMat(3,1) = 2*b*d + 2*a*c ;
  RotMat(3,2) = 2*c*d - 2*a*b ;
  RotMat(3,3) = a^2 - b^2 - c^2 + d^2 ;