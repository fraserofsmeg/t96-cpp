#include <gtest/gtest.h>
extern "C" {
    void t96_(int *iopt, double *parmod, double *ps,
              double *x, double *y, double *z,
              double *bx, double *by, double *bz);
}

void t96(int iopt, const std::array<double,10>& parmod, double ps,
         double x, double y, double z,
         double& bx, double& by, double& bz);

TEST(T96Compare, RandomSamples) {
    const int numTests = 50;
    double parmod[10] = {0};
    for (int i=0;i<numTests;++i) {
        double x = -20.0 + (rand() / (double)RAND_MAX)*30.0;
        double y = -10.0 + (rand() / (double)RAND_MAX)*20.0;
        double z = -10.0 + (rand() / (double)RAND_MAX)*20.0;
        double psi = (-35.0 + (rand() / (double)RAND_MAX)*70.0) * M_PI/180.0;
        int iopt = 1 + rand()%7;
        double bx_f, by_f, bz_f;
        double bx_c, by_c, bz_c;
        t96_(&iopt, parmod, &psi, &x, &y, &z, &bx_f, &by_f, &bz_f);
        std::array<double,10> p = {};
        t96(iopt, p, psi, x, y, z, bx_c, by_c, bz_c);
        EXPECT_NEAR(bx_f, bx_c, 1e-5);
        EXPECT_NEAR(by_f, by_c, 1e-5);
        EXPECT_NEAR(bz_f, bz_c, 1e-5);
    }
}

int main(int argc,char** argv){
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
