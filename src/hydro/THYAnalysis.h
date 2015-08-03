#include<iostream>
#include<fstream>

class THydro3DBj ;
class TRNavier3DBj ;
class TGrid3D ;
class TIC3D ;

class THYAnalysis{

   public:

      virtual ~THYAnalysis() { ; }
      //! Called at the start of a particular hydro run, to perhaps open new
      //! data files or something like that, though often the constructor
      //! handles this job.
      virtual void start() =0;
      //! Called at the end of a particular hydro run
      virtual void stop() =0;
      //! Analyze at time t. The analysis code may decide to print noathing
      //! at this time. This is called each time step.
      virtual void analyze() =0 ;
      //! Called after all runs are completed. Normally this is not needed. If
      //! the destructor can handle the job.
      virtual void terminate() {; }
} ;


//! Purpose:
//!
//! To write the data to make a contour plot of energy density  in the x vs. y
//! or x vs. z or y vs. z plane
//!
//! Can also make similar plots for energy temperature, flow velocity,
//! $\gamma$, and entropy with the output file.
//!
//! Typical usage 
//!
//! double zslice = 1.0;
//! THYAnalysis_2dslice *hya = new  THYAnalysis_2dslice(&rn, &hy, tout, "xy",
//!                                                     zslice) ;
//! 
//! Calling analyze() will produce a data file such as
//!
//! nametag_2dslice_xy_0010.out
//!
//! Which produces slices of the xy plane at a specified z = zslice.
//!
//! The collumns of the output file are
//!
//! t, x, y, z, e, n, p, s, temper, mu, u0, u1, u2, u3, pi11, pi12, pi13, pi22,
//! pi23, pi33 ;
class THYAnalysis_2dslice : public THYAnalysis {

   private:

      double fTNextPrint; //! Next time when a call to analyze() will print.
      double fTOut;       //! Time between printings
      std::string fPlane; //! One of xy, xz, yz
      int fLocation;      //! The index of the third coordinates

      THydro3DBj *fHydro;  
      bool fPlotSolution ; 
      std::ofstream fOut ;

      void write_2dslice()  ;

   public :

      //! Initialize the analysis for a specified plane. 
      //! For definiteness, we take the xy plane. 
      //!
      //! Inputs:
      //!
      //! - time_out Time between printing
      //!
      //! - plane    Specify the orientation of the slize, "xy", "xz", "yz"
      //!
      //! - slice    Specify the third coordinate, e.g. if plane "xy" this
      //!            specifies the z-value of the slice.
      //!
      //! - plot_solution   If plot_solution is true then the initial condition
      //!                   object is called to plot the solution along side
      //!                   the grid.
      //!             
      THYAnalysis_2dslice(THydro3DBj *hy, 
            const double &time_out, 
            const std::string &plane, 
            const double &zslice, 
            const bool &plot_solution=false)  ;

      ~THYAnalysis_2dslice() {; }
      virtual void start() override ;
      virtual void stop() override;
      virtual void analyze() override ;
} ;

//! Purpose:
//!
//! To write a line plot along a coordinate axis to an output appropriate
//! for gnuplot.
//!
//! Can also make similar plots for energy temperature, flow velocity,
//! $\gamma$, and entropy with the output file.
//!
//! Typical usage 
//!
//! double yslice = 4.0;
//! double zslice = 1.0;
//! THYAnalysis_1dslice *hya = new  THYAnalysis_1dslice(&rn, &hy, tout, "x",
//!                                                     yslice, zslice) ;
//! 
//! Calling analyze() will produce a data file such as
//!
//! nametag_2dslice_x_0040_0010.out
//!
//! Which produces slices of the x line at a specified y = yslice and z = zslice
//!
//! The collumns of the output file are
//!
//! t, x, y, z, e, n, p, s, temper, mu, u0, u1, u2, u3, pi11, pi12, pi13, pi22,
//! pi23, pi33 ;
class THYAnalysis_1dslice : public THYAnalysis {

   private:

      double fTNextPrint; //! Next time when a call to analyze() will print.
      double fTOut;       //! Time between printings
      std::string fLine;  //! One  x, y, z
      int fLocation2;     //! The index of the second coordinate
      int fLocation3;     //! The index of the third coordinates

      THydro3DBj *fHydro;  
      bool fPlotSolution ;
      std::ofstream fOut ;

      void write_1dslice()  ;

   public :

      //! Initialize the analysis for a specified plane. 
      //! For definiteness, we take the xy plane. 
      //!
      //! Inputs:
      //!
      //! - time_out Time between printing
      //!
      //! - line     Specify the orientation of the slize, "x", "y", "z"
      //!
      //! - yslice   Specify the second coordinate, e.g. if plane "x" this
      //!            specifies the y-value of the slice
      //!
      //! - zslice   Specify the third coordinate, e.g. if plane "x" this
      //!            specifies the y-value of the slice
      //!
      //! - plot_solution   If plot_solution is true then the initial condition
      //!                   object is called to plot the solution along side
      //!                   the grid.
      THYAnalysis_1dslice(THydro3DBj *hy, double time_out, const std::string &line, const double &yslice, const double &zslice, const bool &plot_solution=false)  ;

      ~THYAnalysis_1dslice() {; }
      virtual void start() override ;
      virtual void stop() override;
      virtual void analyze() override ;
} ;

class THYAnalysis_point : public THYAnalysis {

   private:

      int fLocation1;     //! The index of the first coordinate
      int fLocation2;     //! The index of the second coordinate
      int fLocation3;     //! The index of the third coordinates

      THydro3DBj *fHydro ;  
      bool fPlotSolution ;
      std::ofstream fOut ;

   public :

      //! Initialize the analysis for a specified plane. 
      //! For definiteness, we take the xy plane. 
      //!
      //! Inputs:
      //!
      //! - x,y,z
      //!
      //! - plot_solution   If plot_solution is true then the initial condition
      //!                   object is called to plot the solution along side
      //!                   the grid.
      THYAnalysis_point(THydro3DBj *hy, 
            const double &x, 
            const double &y, 
            const double &z, 
            const bool &plot_solution=false)  ;

      ~THYAnalysis_point() {; }
      virtual void start() override ;
      virtual void stop() override;
      virtual void analyze() override ;
} ;

//! Plots the (integrated) entropy in the transverse plane for a specific eta
//! slice vs. time. The output is a data file 
//!
//! tag_entropy_vs_time_z0000.out with two columns
//!
//! time  \f$\int t dx\,dy s(t, x, y, z)\f$
//!
//! Where \f$t\f$ is short for \f$\tau\f$ and \f$z\f$ is short for \f$\eta\f$
class THYAnalysis_entropy_vs_time : public THYAnalysis {

   private:

      int fLocation3;  //! The index of the third coordinates
      THydro3DBj *fHydro;  
      std::ofstream fOut;
      void write_entropy_vs_time(const double &tau);

   public :

      //! Initialize the analysis model for a specific zslice.
      THYAnalysis_entropy_vs_time(THydro3DBj *hy, const double &zslice); 
      ~THYAnalysis_entropy_vs_time() {; }
      virtual void start() override ;
      virtual void stop() override;
      virtual void analyze() override ;
} ;
