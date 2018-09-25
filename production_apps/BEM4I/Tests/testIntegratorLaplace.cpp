#include "catch.hpp"

#include "../Settings.h"
#include "../FullMatrix.h"
#include "../SurfaceMesh3D.h"
#include "../BESpace.h"
#include "../BEBilinearFormLaplace1Layer.h"
#include "../BEBilinearFormLaplace2Layer.h"
#include "../BEBilinearFormLaplaceHypersingular.h"

using namespace bem4i;

TEST_CASE( "Element system matrices are computed", "[element]" ) {

	typedef double SC;
	typedef double SCVT;
	typedef int LO;

	SurfaceMesh3D< LO, SC > mesh;
	string filename = "../input/lshape_28.txt";
	mesh.load( filename );
	int order = 4;
	int quadOrder[4] = { order, order, order, order };

	std::cout.precision(15);
	std::cout << std::scientific;

	SECTION( "Testing p0p0 integration" ) {

		SC result;
		FullMatrix< LO, SC > elemMatrix( 1, 1 );	
		BESpace< LO, SC > bespaceV( &mesh, p0, p0 );

		SECTION( "Using semi-analytic quadrature" ) {

			quadratureType quadType = Steinbach;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceV, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {
				integrator.computeElemMatrix1Layer( 0, 0, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 8.050529697022868e-02 ).epsilon( 1e-14 ) );
			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrix1Layer( 0, 1, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.830554878450952e-02 ).epsilon( 1e-14 ) ); 

				integrator.computeElemMatrix1Layer( 0, 4, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.711608453472868e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 19, 25, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.711608453472869e-02 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrix1Layer( 0, 2, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.103342514152386e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 0, 5, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.573278126528793e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 6, 17, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.573278126528794e-02 ).epsilon( 1e-14 ) );

			}

		}

		SECTION( "Using fully numerical quadrature" ) {

			quadratureType quadType = SauterSchwab;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceV, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {
				integrator.computeElemMatrix1Layer( 0, 0, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 7.980853550151068e-02 ).epsilon( 1e-14 ) );
			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrix1Layer( 0, 1, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.847765191055143e-02 ).epsilon( 1e-14 ) ); 

				integrator.computeElemMatrix1Layer( 0, 4, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.713972124199284e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 19, 25, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 3.714425417395063e-02 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrix1Layer( 0, 2, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.103408305546696e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 0, 5, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.572450109141559e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix1Layer( 6, 17, elemMatrix );
				result = elemMatrix.get( 0, 0 );
				REQUIRE( result == Approx( 2.572450109141558e-02 ).epsilon( 1e-14 ) );

			}

		}

	}

	SECTION( "Testing p1p0 integration" ) {

		FullMatrix< LO, SC > elemMatrix( 1, 3 );	
		BESpace< LO, SC > bespaceK( &mesh, p1, p0 );

		SECTION( "Using semi-analytic quadrature" ) {

			quadratureType quadType = Steinbach;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceK, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {

				integrator.computeElemMatrix2Layer( 0, 0, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrix2Layer( 0, 1, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 0, 4, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.672476730025871e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.096898356255206e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.648993494407460e-03 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 19, 25, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.096898356255205e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.672476730025870e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.648993494407457e-03 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrix2Layer( 0, 2, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 0, 5, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -8.852093280961244e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.480533852125474e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -3.479614027251039e-03 ).epsilon( 1e-14 ) );


				integrator.computeElemMatrix2Layer( 6, 17, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -8.852093280961249e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.480533852125474e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -3.479614027251036e-03 ).epsilon( 1e-14 ) );

			}

		}

		SECTION( "Using fully numerical quadrature" ) {

			quadratureType quadType = SauterSchwab;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceK, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {

				integrator.computeElemMatrix2Layer( 0, 0, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrix2Layer( 0, 1, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 0, 4, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.692964852243636e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.128733992703048e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.382469235257322e-03 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 19, 25, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.132490009262535e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.700037909588355e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.387198324567223e-03 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrix2Layer( 0, 2, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrix2Layer( 0, 5, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -8.570194662118657e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.489170270776185e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -3.485406471287686e-03 ).epsilon( 1e-14 ) );


				integrator.computeElemMatrix2Layer( 6, 17, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -8.570194662118658e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.489170270776185e-03 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -3.485406471287686e-03 ).epsilon( 1e-14 ) );

			}

		}

	}

	SECTION( "Testing p1p1 integration" ) {

		FullMatrix< LO, SC > elemMatrix( 3, 3 );	
		BESpace< LO, SC > bespaceD( &mesh, p1, p1 );

		SECTION( "Using semi-analytic quadrature" ) {

			quadratureType quadType = Steinbach;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceD, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {

				integrator.computeElemMatrixHypersingular( 0, 0, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 1.610105939404574e-01 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -8.050529697022868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -8.050529697022868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -8.050529697022868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 8.050529697022868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -8.050529697022868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 8.050529697022868e-02 ).epsilon( 1e-14 ) );


			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrixHypersingular( 0, 1, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 3.830554878450952e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 3.830554878450952e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.661109756901904e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -3.830554878450952e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 3.830554878450952e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( -3.830554878450952e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 3.830554878450952e-02 ).epsilon( 1e-14 ) );



				integrator.computeElemMatrixHypersingular( 0, 4, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrixHypersingular( 19, 25, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -3.711608453472868e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 3.711608453472868e-02 ).epsilon( 1e-14 ) );


			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrixHypersingular( 0, 2, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 4.206685028304773e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.103342514152387e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -2.103342514152387e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -2.103342514152387e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 2.103342514152387e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -2.103342514152387e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 2.103342514152387e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrixHypersingular( 0, 5, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.573278126528793e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 2.573278126528793e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 2.573278126528793e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -2.573278126528793e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );


				integrator.computeElemMatrixHypersingular( 6, 17, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 2.573278126528794e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -2.573278126528794e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -2.573278126528794e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( -0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 2.573278126528794e-02 ).epsilon( 1e-14 ) );

			}

		}

		SECTION( "Using fully numerical quadrature" ) {

			quadratureType quadType = SauterSchwab;
			BEIntegratorLaplace< LO, SC > integrator( 
					&bespaceD, quadOrder, quadType);

			SECTION( "Identical triangular panels" ) {

				integrator.computeElemMatrixHypersingular( 0, 0, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 1.596170710030214e-01 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -7.980853550151068e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.980853550151068e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -7.980853550151068e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 7.980853550151068e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -7.980853550151068e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 7.980853550151068e-02 ).epsilon( 1e-14 ) );

			}

			SECTION( "Triangular panels with common edge" ) {

				integrator.computeElemMatrixHypersingular( 0, 1, elemMatrix );

				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 3.847765191055143e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 3.847765191055143e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -7.695530382110287e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -3.847765191055143e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 3.847765191055143e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( -3.847765191055143e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 3.847765191055143e-02 ).epsilon( 1e-14 ) );



				integrator.computeElemMatrixHypersingular( 0, 4, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -3.713972124199284e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 3.713972124199284e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 3.713972124199284e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -3.713972124199284e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrixHypersingular( 19, 25, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 3.714425417395063e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -3.714425417395063e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -3.714425417395063e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 3.714425417395063e-02 ).epsilon( 1e-14 ) );


			}

			SECTION( "Triangular panels with common vertex" ) {

				integrator.computeElemMatrixHypersingular( 0, 2, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 4.206816611093393e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -2.103408305546696e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -2.103408305546696e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( -2.103408305546696e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 2.103408305546696e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -2.103408305546696e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 2.103408305546696e-02 ).epsilon( 1e-14 ) );

				integrator.computeElemMatrixHypersingular( 0, 5, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( -2.572450109141559e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( -0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( 2.572450109141559e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 2.572450109141559e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( -2.572450109141559e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );


				integrator.computeElemMatrixHypersingular( 6, 17, elemMatrix );
				REQUIRE( elemMatrix.get( 0, 0 ) == 
						Approx( 2.572450109141558e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 0, 2 ) == 
						Approx( -2.572450109141558e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 0 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 1 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 1, 2 ) == 
						Approx( 0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 0 ) == 
						Approx( -2.572450109141558e-02 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 1 ) == 
						Approx( -0.000000000000000e+00 ).epsilon( 1e-14 ) );
				REQUIRE( elemMatrix.get( 2, 2 ) == 
						Approx( 2.572450109141558e-02 ).epsilon( 1e-14 ) );



			}

		}

	}

}
