#ifndef __MATRIX_H__
#define __MATRIX_H__


#include<vector>
#include<cassert>
#include<iostream>
namespace Matrix {    
    typedef std::vector<double> row;
    typedef std::vector<row> matrix;

    std::ostream& operator<<(std::ostream& os, const matrix &m) {
        for (int i=0; i<(int)m.size(); i++) {
            for (int j=0; j<(int)m[i].size(); j++) {
                os<<m[i][j];
                if (j+1 != (int)m[i].size()) os<<", ";
            }
            os<<std::endl;
        }
        return os;
    }

    ///Calculate a*b
    matrix operator*(const matrix &a, const matrix &b) {
        assert((int)a.size() > 0 && (int)b.size() > 0 && 
               (int)b[0].size() > 0 && a[0].size() == b.size());
        int n = (int)a.size();
        int p = (int)b.size();
        int m = (int)b[0].size();

        matrix ans(n, row(m));
        for (int i=0; i<n; i++)
            for (int j=0; j<m; j++)
                for (int k=0; k<p; k++)
                    ans[i][j] = (ans[i][j] + a[i][k]*b[k][j]);
        return ans;
    }

    ///Calculate a*b
    matrix operator*(const matrix &a, double d) {
        matrix ans = a;
        for (int i=0; i<(int)ans.size(); i++)
            for (int j=0; j<(int)ans[i].size(); j++)
                 ans[i][j] *= d;
        return ans;
    }

    ///Construct n by n unit matrix
    matrix unitMatrix(int n) {
        matrix ans(n, row(n));
        for (int i=0; i<n; i++)
            ans[i][i] = 1;
        return ans;
    }
    
    ///Construct n by n zero matrix
    matrix zeroMatrix(int n, int m) {
        matrix ans(n, row(m, 0));
        return ans;
    }

    ///calculate a^n
    matrix power(const matrix &a, long long p) {
        if (p == 0)     return unitMatrix((int)a.size());
        matrix ans = power(a, p/2);
        ans = ans * ans;
        if (p%2)        ans = ans*a;
        return ans;
    }

} // namespace matrix



#endif // __MATRIX_H__