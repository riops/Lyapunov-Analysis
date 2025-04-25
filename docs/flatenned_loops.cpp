std::vector<long double> averagedEquationsPolarizationBasisParallel(const std::vector<long double>& allVectors){
  long double mu1 = 1.0;
  long double mu2 = 10.0;

  int N = 2;
  int Imax = 4 * (N * N); // Number of averaged variables
  int CovMax = Imax * Imax; // Number of covariance terms

  // Extract mean vectors
  std::vector<long double> X(allVectors.begin(), allVectors.begin() + Imax);
  std::vector<long double> P(allVectors.begin() + Imax, allVectors.begin() + 2 * Imax);

  // Extract covariance vectors (centered)
  std::vector<long double> XX(allVectors.begin() + 2 * Imax, allVectors.begin() + 2 * Imax + CovMax);
  std::vector<long double> PP(allVectors.begin() + 2 * Imax + CovMax, allVectors.begin() + 2 * Imax + 2 * CovMax);
  std::vector<long double> XP(allVectors.begin() + 2 * Imax + 2 * CovMax, allVectors.begin() + 2 * Imax + 3 * CovMax);

  std::vector<long double> result(2 * Imax + 3 * CovMax);

  // Equations for mean X and P (as before)
  for (int idx_ai_l1m1 = 0; idx_ai_l1m1 < 4*N*N; idx_ai_l1m1++){
    result[idx_ai_l1m1] = P[idx_ai_l1m1];
    for (int lm3 = 0; lm3 < N*N; lm3++){
      for (int lm4 = 0; lm4 < N*N; lm4++){
        for (int lm5 = 0; lm5 < N*N; lm5++){
          for (int idx_kbcd = 1; idx_kbcd <= 16; idx_kbcd++){
            int temp1 = idx_ai_l1m1 % N*N;
            int l1 = static_cast<int>(sqrt(temp1));
            int m1 = temp1 - (l1*l1 + l1);
            int idx_ai = idx_ai_l1m1 / N*N;
            int a = idx_ai / 2 + 1;
            int i = idx_ai % 2 + 1;

            int k = (idx_kbcd - 1) % 2 + 1;
            int b = ((idx_kbcd - 1) / 2) % 2 + 1;
            int c = ((idx_kbcd - 1) / 4) % 2 + 1;
            int d = ((idx_kbcd - 1) / 8) % 2 + 1;

            int l3 = static_cast<int>(sqrt(lm3));
            int m3 = lm3 - (l3*l3 + l3);
            int l4 = static_cast<int>(sqrt(lm4));
            int m4 = lm4 - (l4*l4 + l4);
            int l5 = static_cast<int>(sqrt(lm5));
            int m5 = lm5 - (l5*l5 + l5);

            int idx_bk_l5m5 = indexXInt(b, k, l5, m5, N);
            int idx_ck_di_l3m3_l4m4 = indexXXInt(c, k, l3, m3, d, i, l4, m4, N);
            int idx_ck_l3m3 = indexXInt(c, k, l3, m3, N);
            int idx_bk_di_l5m5_l4m4 = indexXXInt(b, k, l5, m5, d, i, l4, m4, N);
            int idx_di_l4m4 = indexXInt(d, i, l4, m4, N);
            int idx_ck_bk_l3m3_l5m5 = indexXXInt(c, k, l3, m3, b, k, l5, m5, N);

            long double pDotSubSum = HFunctionInt(l3, l4, l5, l1, m3, m4, m5, m1, N) * GInt(a, b, c, d) * ( X[idx_bk_l5m5]*XX[idx_ck_di_l3m3_l4m4] + X[idx_ck_l3m3]*XX[idx_bk_di_l5m5_l4m4] + X[idx_di_l4m4]*XX[idx_ck_bk_l3m3_l5m5] + X[idx_bk_l5m5] * X[idx_ck_l3m3] * X[idx_di_l4m4]);
            result[Imax + idx_ai_l1m1] = -((mu1 + l1*(l1 + 1)) * X[idx_ai_l1m1] - mu2 * pDotSubSum);

            for (int idx_ej_l2m2 = 0; idx_ej_l2m2 < 4*N*N; idx_ej_l2m2++){
              int temp2 = idx_ej_l2m2 % N*N;
              int l2 = static_cast<int>(sqrt(temp2));
              int m2 = temp2 - (l2*l2 + l2);
              int idx_ej = idx_ej_l2m2 / N*N;
              int e = idx_ej / 2 + 1;
              int j = idx_ej % 2 + 1;

              int idx_ai_ej_l1m1_l2m2 = indexXXInt(a, i, l1, m1, e, j, l2, m2, N);
              // First term of the XXdot term
              result[2 * Imax + idx_ai_ej_l1m1_l2m2] = XP[idx_ai_ej_l1m1_l2m2];
              int idx_bk_ej_l5m5_l2m2 = indexXXInt(b, k, l5, m5, e, j, l2, m2, N);
              int idx_ck_ej_l3m3_l2m2 = indexXXInt(c, k, l3, m3, e, j, l2, m2, N);
              int idx_di_ej_l4m4_l2m2 = indexXXInt(d, i, l4, m4, e, j, l2, m2, N);
              // {i, a, l1, m1} <-> {j, e, l2, m2}
              int idx_ej_ai_l2m2_l1m1 = indexXXInt(e, j, l2, m2, a, i, l1, m1, N);
              int idx_bk_ai_l5m5_l1m1 = indexXXInt(b, k, l5, m5, a, i, l1, m1, N);
              int idx_ck_dj_l3m3_l4m4 = indexXXInt(c, k, l3, m3, d, j, l4, m4, N);
              int idx_ck_ai_l3m3_l1m1 = indexXXInt(c, k, l3, m3, a, i, l1, m1, N);
              int idx_bk_dj_l5m5_l4m4 = indexXXInt(b, k, l5, m5, d, j, l4, m4, N);
              int idx_dj_ai_l4m4_l1m1 = indexXXInt(d, j, l4, m4, a, i, l1, m1, N);
              int idx_dj_l4m4 = indexXInt(d, j, l4, m4, N);
              // Compute the dot product of the PP term.
              long double ppDotSubSum1 = HFunctionInt(l3, l4, l5, l1, m3, m4, m5, m1, N) * GInt(a, b, c, d) * ( XP[idx_bk_ej_l5m5_l2m2]*XX[idx_ck_di_l3m3_l4m4] + XP[idx_ck_ej_l3m3_l2m2]*XX[idx_bk_di_l5m5_l4m4] + XP[idx_di_ej_l4m4_l2m2]*XX[idx_ck_bk_l3m3_l5m5] + XP[idx_bk_ej_l5m5_l2m2]*X[idx_ck_l3m3]*X[idx_di_l4m4] + XP[idx_ck_ej_l3m3_l2m2]*X[idx_bk_l5m5]*X[idx_di_l4m4] + XP[idx_di_ej_l4m4_l2m2]*X[idx_ck_l3m3]*X[idx_bk_l5m5] );
              long double ppDotSubSum2 = HFunctionInt(l3, l4, l5, l2, m3, m4, m5, m2, N) * GInt(e, b, c, d) * ( XP[idx_bk_ai_l5m5_l1m1]*XX[idx_ck_dj_l3m3_l4m4] + XP[idx_ck_ai_l3m3_l1m1]*XX[idx_bk_dj_l5m5_l4m4] + XP[idx_dj_ai_l4m4_l1m1]*XX[idx_ck_bk_l3m3_l5m5] + XP[idx_bk_ai_l5m5_l1m1]*X[idx_ck_l3m3]*X[idx_dj_l4m4] + XP[idx_ck_ai_l3m3_l1m1]*X[idx_bk_l5m5]*X[idx_dj_l4m4] + XP[idx_dj_ai_l4m4_l1m1]*X[idx_ck_l3m3]*X[idx_bk_l5m5] );

              result[2 * Imax + CovMax + idx_ai_ej_l1m1_l2m2] = -((mu1 + l1*(l1 + 1) + l2*(l2 + 1)) * XP[idx_ai_ej_l1m1_l2m2] - mu2 * (ppDotSubSum1 + ppDotSubSum2));

              //Compute the dot product of the XP term.
              long double xpDotSubSum = HFunctionInt(l3, l4, l5, l1, m3, m4, m5, m1, N) * GInt(a, b, c, d) * ( XX[idx_bk_ej_l5m5_l2m2]*XX[idx_ck_di_l3m3_l4m4] + XX[idx_ck_ej_l3m3_l2m2]*XX[idx_bk_di_l5m5_l4m4] + XX[idx_di_ej_l4m4_l2m2]*XX[idx_ck_bk_l3m3_l5m5] + XX[idx_bk_ej_l5m5_l2m2]*X[idx_ck_l3m3]*X[idx_di_l4m4] + XX[idx_ck_ej_l3m3_l2m2]*X[idx_bk_l5m5]*X[idx_di_l4m4] + XX[idx_di_ej_l4m4_l2m2]*X[idx_ck_l3m3]*X[idx_bk_l5m5] );

              result[2 * Imax + 2 * CovMax + idx_ai_ej_l1m1_l2m2] = -((mu1 + l1*(l1 + 1)) * XX[idx_ai_ej_l1m1_l2m2] - mu2 * xpDotSubSum) + PP[idx_ai_ej_l1m1_l2m2];
            }
          }
        }
      }
    }
  }
  return result;
}
