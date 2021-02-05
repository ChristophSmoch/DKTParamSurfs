#ifndef __IPOPTINCLUDES_H
#define __IPOPTINCLUDES_H

#ifdef PESOPT_WITH_IPOPT

#ifdef __GNUC__
#pragma GCC system_header
#endif

#include <cstddef>
#include <pesopt_IO.h>
#include <energyDefines.h>


#define HAVE_CSTDDEF
#include <IpTNLP.hpp>
#undef HAVE_CSTDDEF

#include "IpIpoptApplication.hpp"
#include "IpSolveStatistics.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpOrigIpoptNLP.hpp"



template <typename RealType>
void setToleranceOptionsIpopt ( const RealType ipoptTol, 
                                const int MaxIterations, 
                                Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp ){
    ipoptApp->Options()->SetNumericValue ( "tol", ipoptTol * 1e-2 );
    ipoptApp->Options()->SetNumericValue ( "acceptable_tol", ipoptTol );
    //ipoptApp->Options()->SetNumericValue( "acceptable_constr_viol_tol", 1.e-6 );
    ipoptApp->Options ()->SetIntegerValue ( "max_iter", MaxIterations );
}

void setPrintOptionsIpopt ( const int ipoptPrintLevel, 
                            const string ipoptOutputFile,
                            Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp ){
    ipoptApp->Options()->SetIntegerValue( "print_level", ipoptPrintLevel );
    if ( ipoptOutputFile != "" ){
        ipoptApp->Options()->SetIntegerValue( "file_print_level", ipoptPrintLevel );
        ipoptApp->Options()->SetStringValue( "output_file", ipoptOutputFile.c_str() );
    }   
}

void switchLinearSolverTypeIpopt ( const int linearSolverTypeIpopt,
                                   Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp ) {
    switch( linearSolverTypeIpopt ){
        case  0: ipoptApp->Options()->SetStringValue ( "linear_solver", "MUMPS" ); break;
        case 27: ipoptApp->Options()->SetStringValue ( "linear_solver", "ma27" ); break;
        case 57: ipoptApp->Options()->SetStringValue ( "linear_solver", "ma57" ); break;
        case 77: ipoptApp->Options()->SetStringValue ( "linear_solver", "ma77" ); break;
        case 86: ipoptApp->Options()->SetStringValue ( "linear_solver", "ma86" ); break;
        case 97: ipoptApp->Options()->SetStringValue ( "linear_solver", "ma97" ); break;
        default: break;
    }
}

void setDerivateTestOptionsIpopt ( Ipopt::SmartPtr< Ipopt::IpoptApplication > ipoptApp ){
    ipoptApp->Options()->SetStringValue  ( "derivative_test_print_all", "yes" );
    ipoptApp->Options()->SetNumericValue ( "derivative_test_tol", 1e-4 );
    ipoptApp->Options()->SetNumericValue ( "derivative_test_perturbation", 1e-8 );
    ipoptApp->Options()->SetIntegerValue ( "max_iter", 0 );
}


string getIpoptStatus ( Ipopt::ApplicationReturnStatus ipoptStatus  ) {
    
    string status = "";
    
    switch( ipoptStatus ){
        case 0 :
            status = "Solve Succeeded";
            break;
            
        case 1 :
            status = "Solved To Acceptable Level";
            break;
            
        case 2 :
            status = "Infeasible Problem Detected";
            break;
            
        case 3 :
            status = "Search Direction Becomes Too Small";
            break;
            
        case 4 :
            status = "Diverging Iterates";
            break;
        
        case 5 :
            status = "User Requested Stop";
            break;
            
        case 6 :
            status = "Feasible Point Found";
            break;
            
        case -1 :
            status = "Maximum Iterations Exceeded";
            break;
            
        case -2 :
            status = "Restoration Failed";
            break;
            
        case -3 :
            status = "Error In Step Computation";
            break;
            
        case -10 :
            status = "Not Enough Degrees Of Freedom";
            break;
            
        case -11 :
            status = "Invalid Problem Definition";
            break;
            
        case -12 :
            status = "Invalid Option";
            break;
            
        case -13 :
            status = "Invalid Number Detected";
            break;
            
        case -100 :
            status = "Unrecoverable Exception";
            break;
            
        case -101 :
            status = "NonIpopt Exception Thrown";
            break;
            
        case -102 :
            status = "Insufficient Memory";
            break;
            
        case -199 :
            status = "Internal Error";
            break;
    }
    
    return status;

}


void outputIpoptStatus ( Ipopt::ApplicationReturnStatus ipoptStatus, const bool outputOnlyIfFailed = false ) {
    string status = getIpoptStatus( ipoptStatus );
    // status: 0 - Solve_Succeeded, 1 - Solved_To_Acceptable_Level
    if( !outputOnlyIfFailed && ( ipoptStatus == 0 || ipoptStatus == 1 ) ) cout << "Ipopt finished with status: " << status.c_str() << endl;
}







#endif //PESOPT_WITH_IPOPT

#endif //__IPOPTINCLUDES_H
