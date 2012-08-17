
ThetaPhi := proc (x::(Vector(3))) local theta, phi; theta := arccos(x[3]); if x[1] = 0 and x[2] = 0 then phi := 0 elif x[1] = 0 and 0 < x[2] then phi := (1/2)*evalf(Pi) elif x[1] = 0 and x[2] < 0 then phi := -(1/2)*evalf(Pi) elif x[1] < 0 then phi := arctan(x[2]/x[1])+evalf(Pi) elif 0 < x[1] then phi := arctan(x[2]/x[1]) end if; theta, phi end proc;
%;

GeneralFractionalAnisotropyFunction := proc (B::Matrix, y::(Vector(45)), N::integer)::float; local A, SumOf, SumOfSquares; A := B.y; SumOf := add(A[i], i = 1 .. N); SumOfSquares := add(A[k]^2, k = 1 .. N); return sqrt((N-SumOf^2/SumOfSquares)/(N-1)) end proc;
%;


GeodesicConcentrationFunction := proc (vec::(Vector(45))) local gc; gc := evalf(sqrt(Pi)*(.2*vec[1]+(4/35)*sqrt(5)*vec[4]+(8/315)*vec[11])); return gc end proc;


OdfValue := proc (x::(Vector(3)), y::(Vector(45))) local a, b, l, m; global ans; a := evalf(ThetaPhi(x)); b := [0, 1, 6, 15, 28]; ans := 0; for l from 0 by 2 to 8 do for m from -l to l do if m < 0 then ans := ans+y[l+m+1+b[(1/2)*l+1]]*evalf(sqrt(2))*Re(SphericalY(l, m, a[1], a[2])) end if; if m = 0 then ans := ans+y[l+m+1+b[(1/2)*l+1]]*Re(SphericalY(l, 0, a[1], a[2])) end if; if 0 < m then ans := ans+y[l+m+1+b[(1/2)*l+1]]*evalf(sqrt(2))*Im(SphericalY(l, m, a[1], a[2])) end if end do end do; return ans end proc;

%;
