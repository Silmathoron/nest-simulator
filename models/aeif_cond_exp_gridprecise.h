/*
 *  aeif_cond_exp_gridprecise.h
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef AEIF_COND_EXP_GP_H
#define AEIF_COND_EXP_GP_H

// Generated includes:
#include "config.h"

#ifdef HAVE_GSL

// External includes:
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

// Includes from nestkernel:
#include "archiving_node.h"
#include "connection.h"
#include "event.h"
#include "nest_types.h"
#include "recordables_map.h"
#include "ring_buffer.h"
#include "universal_data_logger.h"

/* BeginDocumentation
Name: aeif_cond_exp_gridprecise - Conductance based exponential integrate-and-
  fire neuron model according to Brette and Gerstner (2005), implementing a
  linear interpolation to find the precise time where the threshold was
  crossed, i.e. the spiking time.

Description:
aeif_cond_alpha is the adaptive exponential integrate and fire neuron according
to Brette and Gerstner (2005) and synaptic conductances are modelled as alpha
functions. This model implements a linear interpolation to find spike times
more precisely.

This implementation uses the embedded 4th order Runge-Kutta-Fehlberg solver
with adaptive stepsize to integrate the differential equation.

The membrane potential is given by the following differential equation:

C dV/dt = -g_L*(V-E_L) + g_L*Delta_T*exp((V-V_T)/Delta_T) - g_e(t)*(V-E_e)
          -g_i(t)*(V-E_i) - w + I_e

and

tau_w * dw/dt = a*(V-E_L) - w

Parameters:
The following parameters can be set in the status dictionary.

Dynamic state variables:
  V_m        double - Membrane potential in mV
  g_ex       double - Excitatory synaptic conductance in nS.
  dg_ex      double - First derivative of g_ex in nS/ms
  g_in       double - Inhibitory synaptic conductance in nS.
  dg_in      double - First derivative of g_in in nS/ms.
  w          double - Spike-adaptation current in pA.

Membrane Parameters:
  C_m        double - Capacity of the membrane in pF
  t_ref      double - Duration of refractory period in ms.
  V_reset    double - Reset value for V_m after a spike. In mV.
  E_L        double - Leak reversal potential in mV.
  g_L        double - Leak conductance in nS.
  I_e        double - Constant external input current in pA.

Spike adaptation parameters:
  a          double - Subthreshold adaptation in nS.
  b          double - Spike-triggered adaptation in pA.
  Delta_T    double - Slope factor in mV
  tau_w      double - Adaptation time constant in ms
  V_th       double - Spike initiation threshold in mV
  V_peak     double - Spike detection threshold in mV.

Synaptic parameters
  E_ex       double - Excitatory reversal potential in mV.
  tau_syn_ex double - Characteristic decrease time of excitatory synaptic conductance in ms
(exponential function).
  E_in       double - Inhibitory reversal potential in mV.
  tau_syn_in double - Characteristic decrease time of inhibitory synaptic conductance in ms
(exponential function).

Integration parameters
  gsl_error_tol  double - This parameter controls the admissible error of the GSL integrator.
                          Reduce it if NEST complains about numerical instabilities.

Author: Tanguy Fardet, modified from Marc-Oliver Gewaltig's implementation

Sends: SpikeEvent

Receives: SpikeEvent, CurrentEvent, DataLoggingRequest

References: Brette R and Gerstner W (2005) Adaptive Exponential Integrate-and-
  Fire Model as an Effective Description of Neuronal Activity.
  J Neurophysiol 94:3637-3642

SeeAlso: iaf_cond_alpha, aeif_cond_exp, aeif_cond_alpha_gridprecise,
*/

namespace nest
{
/**
 * Function computing right-hand side of ODE for GSL solver.
 * @note Must be declared here so we can befriend it in class.
 * @note Must have C-linkage for passing to GSL. Internally, it is
 *       a first-class C++ function, but cannot be a member function
 *       because of the C-linkage.
 * @note No point in declaring it inline, since it is called
 *       through a function pointer.
 * @param void* Pointer to model neuron instance.
 */
extern "C" int aeif_cond_exp_gridprecise_dynamics( double, const double*, double*, void* );

class aeif_cond_exp_gridprecise : public Archiving_Node
{

public:
  aeif_cond_exp_gridprecise();
  aeif_cond_exp_gridprecise( const aeif_cond_exp_gridprecise& );
  ~aeif_cond_exp_gridprecise();

  /**
   * Import sets of overloaded virtual functions.
   * @see Technical Issues / Virtual Functions: Overriding, Overloading, and Hiding
   */
  using Node::handle;
  using Node::handles_test_event;

  port send_test_event( Node&, rport, synindex, bool );

  void handle( SpikeEvent& );
  void handle( CurrentEvent& );
  void handle( DataLoggingRequest& );

  port handles_test_event( SpikeEvent&, rport );
  port handles_test_event( CurrentEvent&, rport );
  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

private:
  void init_state_( const Node& proto );
  void init_buffers_();
  void calibrate();
  void update( Time const&, const long_t, const long_t );
  void interpolate_( double&, double );
  void spiking_( Time const&, const long_t, const double );

  // END Boilerplate function declarations ----------------------------

  // Friends --------------------------------------------------------

  // make dynamics function quasi-member
  friend int aeif_cond_exp_gridprecise_dynamics( double, const double*, double*, void* );

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< aeif_cond_exp_gridprecise >;
  friend class UniversalDataLogger< aeif_cond_exp_gridprecise >;

private:
  // ----------------------------------------------------------------

  //! Independent parameters
  struct Parameters_
  {
    double_t V_peak_;  //!< Spike detection threshold in mV
    double_t V_reset_; //!< Reset Potential in mV
    double_t t_ref_;   //!< Refractory period in ms

    double_t g_L;        //!< Leak Conductance in nS
    double_t C_m;        //!< Membrane Capacitance in pF
    double_t E_ex;       //!< Excitatory reversal Potential in mV
    double_t E_in;       //!< Inhibitory reversal Potential in mV
    double_t E_L;        //!< Leak reversal Potential (aka resting potential) in mV
    double_t Delta_T;    //!< Slope faktor in ms.
    double_t tau_w;      //!< adaptation time-constant in ms.
    double_t a;          //!< Subthreshold adaptation in nS.
    double_t b;          //!< Spike-triggered adaptation in pA
    double_t V_th;       //!< Spike threshold in mV.
    double_t t_ref;      //!< Refractory period in ms.
    double_t tau_syn_ex; //!< Excitatory synaptic rise time.
    double_t tau_syn_in; //!< Excitatory synaptic rise time.
    double_t I_e;        //!< Intrinsic current in pA.

    double_t gsl_error_tol; //!< error bound for GSL integrator

    Parameters_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    void set( const DictionaryDatum& ); //!< Set values from dicitonary
  };

public:
  // ----------------------------------------------------------------

  /**
   * State variables of the model.
   * @note Copy constructor and assignment operator required because
   *       of C-style array.
   */
  struct State_
  {
    /**
     * Enumeration identifying elements in state array State_::y_.
     * The state vector must be passed to GSL as a C array. This enum
     * identifies the elements of the vector. It must be public to be
     * accessible from the iteration function.
     */
    enum StateVecElems
    {
      V_M = 0,
      G_EXC, // 1
      G_INH, // 2
      W,     // 3
      STATE_VEC_SIZE
    };

    double_t y_[ STATE_VEC_SIZE ];     //!< neuron state, must be C-array for GSL solver
    double_t y_old_[ STATE_VEC_SIZE ]; //!< old neuron state, must be C-array for GSL solver
    int_t r_;                          //!< number of refractory steps remaining
    double_t r_offset_; // offset on the refractory time if it is not a multiple of step_

    State_( const Parameters_& ); //!< Default initialization
    State_( const State_& );
    State_& operator=( const State_& );

    void get( DictionaryDatum& ) const;
    void set( const DictionaryDatum&, const Parameters_& );
  };

  // ----------------------------------------------------------------

  /**
   * Buffers of the model.
   */
  struct Buffers_
  {
    Buffers_( aeif_cond_exp_gridprecise& );                  //!<Sets buffer pointers to 0
    Buffers_( const Buffers_&, aeif_cond_exp_gridprecise& ); //!<Sets buffer pointers to 0

    //! Logger for all analog data
    UniversalDataLogger< aeif_cond_exp_gridprecise > logger_;

    /** buffers and sums up incoming spikes/currents */
    RingBuffer spike_exc_;
    RingBuffer spike_inh_;
    RingBuffer currents_;

    /** GSL ODE stuff */
    gsl_odeiv_step* s_;    //!< stepping function
    gsl_odeiv_control* c_; //!< adaptive stepsize control function
    gsl_odeiv_evolve* e_;  //!< evolution function
    gsl_odeiv_system sys_; //!< struct describing system

    // IntergrationStep_ should be reset with the neuron on ResetNetwork,
    // but remain unchanged during calibration. Since it is initialized with
    // step_, and the resolution cannot change after nodes have been created,
    // it is safe to place both here.
    double_t step_;          //!< step size in ms
    double IntegrationStep_; //!< current integration time step, updated by GSL

    /**
     * Input current injected by CurrentEvent.
     * This variable is used to transport the current applied into the
     * _dynamics function computing the derivative of the state vector.
     * It must be a part of Buffers_, since it is initialized once before
     * the first simulation, but not modified before later Simulate calls.
     */
    double_t I_stim_;
  };

  // ----------------------------------------------------------------

  /**
   * Internal variables of the model.
   */
  struct Variables_
  {
    int_t RefractoryCounts_;
    double_t RefractoryOffset_;
  };

  // Access functions for UniversalDataLogger -------------------------------

  //! Read out state vector elements, used by UniversalDataLogger
  template < State_::StateVecElems elem >
  double_t
  get_y_elem_() const
  {
    return S_.y_[ elem ];
  }

  // ----------------------------------------------------------------

  Parameters_ P_;
  State_ S_;
  Variables_ V_;
  Buffers_ B_;

  //! Mapping of recordables names to access functions
  static RecordablesMap< aeif_cond_exp_gridprecise > recordablesMap_;
};

inline port
aeif_cond_exp_gridprecise::send_test_event( Node& target, rport receptor_type, synindex, bool )
{
  SpikeEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
aeif_cond_exp_gridprecise::handles_test_event( SpikeEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return 0;
}

inline port
aeif_cond_exp_gridprecise::handles_test_event( CurrentEvent&, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return 0;
}

inline port
aeif_cond_exp_gridprecise::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
    throw UnknownReceptorType( receptor_type, get_name() );
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
aeif_cond_exp_gridprecise::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  S_.get( d );
  Archiving_Node::get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
aeif_cond_exp_gridprecise::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d );         // throws if BadProperty
  State_ stmp = S_;      // temporary copy in case of errors
  stmp.set( d, ptmp );   // throws if BadProperty

  // We now know that (ptmp, stmp) are consistent. We do not
  // write them back to (P_, S_) before we are also sure that
  // the properties to be set in the parent class are internally
  // consistent.
  Archiving_Node::set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
  S_ = stmp;
}

} // namespace

#endif // HAVE_GSL
#endif // AEIF_COND_EXP_GP_H
