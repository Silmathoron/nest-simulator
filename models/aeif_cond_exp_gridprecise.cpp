/*
 *  aeif_cond_exp_gridprecise.cpp
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

#include "aeif_cond_exp_gridprecise.h"
#ifdef HAVE_GSL_1_11

// C++ includes:
#include <cmath>
#include <cstdio>
#include <iomanip>
#include <iostream>
#include <limits>

// Includes from libnestutil:
#include "numerics.h"

// Includes from nestkernel:
#include "exceptions.h"
#include "kernel_manager.h"
#include "nest_names.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

/* ----------------------------------------------------------------
 * Recordables map
 * ---------------------------------------------------------------- */

nest::RecordablesMap< nest::aeif_cond_exp_gridprecise >
  nest::aeif_cond_exp_gridprecise::recordablesMap_;

namespace nest // template specialization must be placed in namespace
{
// Override the create() method with one call to RecordablesMap::insert_()
// for each quantity to be recorded.
template <>
void
RecordablesMap< aeif_cond_exp_gridprecise >::create()
{
  // use standard names whereever you can for consistency!
  insert_(
    names::V_m, &aeif_cond_exp_gridprecise::get_y_elem_< aeif_cond_exp_gridprecise::State_::V_M > );
  insert_( names::g_ex,
    &aeif_cond_exp_gridprecise::get_y_elem_< aeif_cond_exp_gridprecise::State_::G_EXC > );
  insert_( names::g_in,
    &aeif_cond_exp_gridprecise::get_y_elem_< aeif_cond_exp_gridprecise::State_::G_INH > );
  insert_(
    names::w, &aeif_cond_exp_gridprecise::get_y_elem_< aeif_cond_exp_gridprecise::State_::W > );
}
}

extern "C" int
nest::aeif_cond_exp_gridprecise_dynamics( double, const double y[], double f[], void* pnode )
{
  // a shorthand
  typedef nest::aeif_cond_exp_gridprecise::State_ S;

  // get access to node so we can almost work as in a member function
  assert( pnode );
  const nest::aeif_cond_exp_gridprecise& node =
    *( reinterpret_cast< nest::aeif_cond_exp_gridprecise* >( pnode ) );

  // y[] here is---and must be---the state vector supplied by the integrator,
  // not the state vector in the node, node.S_.y[].

  // The following code is verbose for the sake of clarity. We assume that a
  // good compiler will optimize the verbosity away ...

  // shorthand for state variables
  const double_t& V = y[ S::V_M ];
  const double_t& g_ex = y[ S::G_EXC ];
  const double_t& g_in = y[ S::G_INH ];
  const double_t& w = y[ S::W ];

  const double_t I_syn_exc = g_ex * ( V - node.P_.E_ex );
  const double_t I_syn_inh = g_in * ( V - node.P_.E_in );

  // We pre-compute the argument of the exponential
  const double_t exp_arg = ( V - node.P_.V_th ) / node.P_.Delta_T;

  // Upper bound for exponential argument to avoid numerical instabilities
  const double_t MAX_EXP_ARG = 10.;

  // If the argument is too large, we clip it.
  const double_t I_spike = node.P_.Delta_T * std::exp( std::min( exp_arg, MAX_EXP_ARG ) );

  // dv/dt
  f[ S::V_M ] = ( -node.P_.g_L * ( ( V - node.P_.E_L ) - I_spike ) - I_syn_exc - I_syn_inh - w
                  + node.P_.I_e + node.B_.I_stim_ ) / node.P_.C_m;

  f[ S::G_EXC ] = -g_ex / node.P_.tau_syn_ex; // Synaptic Conductance (nS)
  f[ S::G_INH ] = -g_in / node.P_.tau_syn_in; // Synaptic Conductance (nS)

  // Adaptation current w.
  f[ S::W ] = ( node.P_.a * ( V - node.P_.E_L ) - w ) / node.P_.tau_w;

  return GSL_SUCCESS;
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameters and state
 * ---------------------------------------------------------------- */

nest::aeif_cond_exp_gridprecise::Parameters_::Parameters_()
  : V_peak_( 0.0 )    // mV, should not be larger that V_th+10
  , V_reset_( -60.0 ) // mV
  , t_ref_( 0.0 )     // ms
  , g_L( 30.0 )       // nS
  , C_m( 281.0 )      // pF
  , E_ex( 0.0 )       // mV
  , E_in( -85.0 )     // mV
  , E_L( -70.6 )      // mV
  , Delta_T( 2.0 )    // mV
  , tau_w( 144.0 )    // ms
  , a( 4.0 )          // nS
  , b( 80.5 )         // pA
  , V_th( -50.4 )     // mV
  , tau_syn_ex( 0.2 ) // ms
  , tau_syn_in( 2.0 ) // ms
  , I_e( 0.0 )        // pA
  , gsl_error_tol( 1e-6 )
{
}

nest::aeif_cond_exp_gridprecise::State_::State_( const Parameters_& p )
  : r_( 0 )
  , r_offset_( 0. )
{
  y_[ 0 ] = p.E_L;
  for ( size_t i = 1; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = 0;
}

nest::aeif_cond_exp_gridprecise::State_::State_( const State_& s )
  : r_( s.r_ )
  , r_offset_( s.r_offset_ )
{
  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
}

nest::aeif_cond_exp_gridprecise::State_& nest::aeif_cond_exp_gridprecise::State_::operator=(
  const State_& s )
{
  assert( this != &s ); // would be bad logical error in program

  for ( size_t i = 0; i < STATE_VEC_SIZE; ++i )
    y_[ i ] = s.y_[ i ];
  r_ = s.r_;
  r_offset_ = s.r_offset_;
  return *this;
}

/* ----------------------------------------------------------------
 * Parameter and state extractions and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::aeif_cond_exp_gridprecise::Parameters_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::C_m, C_m );
  def< double >( d, names::V_th, V_th );
  def< double >( d, names::t_ref, t_ref_ );
  def< double >( d, names::g_L, g_L );
  def< double >( d, names::E_L, E_L );
  def< double >( d, names::V_reset, V_reset_ );
  def< double >( d, names::E_ex, E_ex );
  def< double >( d, names::E_in, E_in );
  def< double >( d, names::tau_syn_ex, tau_syn_ex );
  def< double >( d, names::tau_syn_in, tau_syn_in );
  def< double >( d, names::a, a );
  def< double >( d, names::b, b );
  def< double >( d, names::Delta_T, Delta_T );
  def< double >( d, names::tau_w, tau_w );
  def< double >( d, names::I_e, I_e );
  def< double >( d, names::V_peak, V_peak_ );
  def< double >( d, names::gsl_error_tol, gsl_error_tol );
}

void
nest::aeif_cond_exp_gridprecise::Parameters_::set( const DictionaryDatum& d )
{
  updateValue< double >( d, names::V_th, V_th );
  updateValue< double >( d, names::V_peak, V_peak_ );
  updateValue< double >( d, names::t_ref, t_ref_ );
  updateValue< double >( d, names::E_L, E_L );
  updateValue< double >( d, names::V_reset, V_reset_ );
  updateValue< double >( d, names::E_ex, E_ex );
  updateValue< double >( d, names::E_in, E_in );

  updateValue< double >( d, names::C_m, C_m );
  updateValue< double >( d, names::g_L, g_L );

  updateValue< double >( d, names::tau_syn_ex, tau_syn_ex );
  updateValue< double >( d, names::tau_syn_in, tau_syn_in );

  updateValue< double >( d, names::a, a );
  updateValue< double >( d, names::b, b );
  updateValue< double >( d, names::Delta_T, Delta_T );
  updateValue< double >( d, names::tau_w, tau_w );

  updateValue< double >( d, names::I_e, I_e );

  updateValue< double >( d, names::gsl_error_tol, gsl_error_tol );

  if ( V_peak_ <= V_th )
    throw BadProperty( "V_peak must be larger than threshold." );

  if ( V_reset_ >= V_peak_ )
    throw BadProperty( "Ensure that: V_reset < V_peak ." );

  if ( C_m <= 0 )
    throw BadProperty( "Capacitance must be strictly positive." );

  if ( g_L <= 0 )
    throw BadProperty( "Leak conductance must be strictly positive." );

  if ( t_ref_ < 0 )
    throw BadProperty( "Refractory time cannot be negative." );

  if ( tau_syn_ex <= 0 || tau_syn_in <= 0 || tau_w <= 0 )
    throw BadProperty( "All time constants must be strictly positive." );

  if ( gsl_error_tol <= 0. )
    throw BadProperty( "The gsl_error_tol must be strictly positive." );
}

void
nest::aeif_cond_exp_gridprecise::State_::get( DictionaryDatum& d ) const
{
  def< double >( d, names::V_m, y_[ V_M ] );
  def< double >( d, names::g_ex, y_[ G_EXC ] );
  def< double >( d, names::g_in, y_[ G_INH ] );
  def< double >( d, names::w, y_[ W ] );
}

void
nest::aeif_cond_exp_gridprecise::State_::set( const DictionaryDatum& d, const Parameters_& )
{
  updateValue< double >( d, names::V_m, y_[ V_M ] );
  updateValue< double >( d, names::g_ex, y_[ G_EXC ] );
  updateValue< double >( d, names::g_in, y_[ G_INH ] );
  updateValue< double >( d, names::w, y_[ W ] );

  if ( y_[ G_EXC ] < 0 || y_[ G_INH ] < 0 )
    throw BadProperty( "Conductances must not be negative." );
}

nest::aeif_cond_exp_gridprecise::Buffers_::Buffers_( aeif_cond_exp_gridprecise& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

nest::aeif_cond_exp_gridprecise::Buffers_::Buffers_( const Buffers_&, aeif_cond_exp_gridprecise& n )
  : logger_( n )
  , s_( 0 )
  , c_( 0 )
  , e_( 0 )
{
  // Initialization of the remaining members is deferred to
  // init_buffers_().
}

/* ----------------------------------------------------------------
 * Default and copy constructor for node, and destructor
 * ---------------------------------------------------------------- */

nest::aeif_cond_exp_gridprecise::aeif_cond_exp_gridprecise()
  : Archiving_Node()
  , P_()
  , S_( P_ )
  , B_( *this )
{
  recordablesMap_.create();
}

nest::aeif_cond_exp_gridprecise::aeif_cond_exp_gridprecise( const aeif_cond_exp_gridprecise& n )
  : Archiving_Node( n )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}

nest::aeif_cond_exp_gridprecise::~aeif_cond_exp_gridprecise()
{
  // GSL structs may not have been allocated, so we need to protect destruction
  if ( B_.s_ )
    gsl_odeiv_step_free( B_.s_ );
  if ( B_.c_ )
    gsl_odeiv_control_free( B_.c_ );
  if ( B_.e_ )
    gsl_odeiv_evolve_free( B_.e_ );
}

/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::aeif_cond_exp_gridprecise::init_state_( const Node& proto )
{
  const aeif_cond_exp_gridprecise& pr = downcast< aeif_cond_exp_gridprecise >( proto );
  S_ = pr.S_;
}

void
nest::aeif_cond_exp_gridprecise::init_buffers_()
{
  B_.spike_exc_.clear(); // includes resize
  B_.spike_inh_.clear(); // includes resize
  B_.currents_.clear();  // includes resize
  Archiving_Node::clear_history();

  B_.logger_.reset();

  B_.step_ = Time::get_resolution().get_ms();

  // We must integrate this model with high-precision to obtain decent results
  B_.IntegrationStep_ = std::min( 0.01, B_.step_ );

  if ( B_.s_ == 0 )
    B_.s_ = gsl_odeiv_step_alloc( gsl_odeiv_step_rkf45, State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_step_reset( B_.s_ );

  if ( B_.c_ == 0 )
    B_.c_ = gsl_odeiv_control_yp_new( P_.gsl_error_tol, P_.gsl_error_tol );
  else
    gsl_odeiv_control_init( B_.c_, P_.gsl_error_tol, P_.gsl_error_tol, 0.0, 1.0 );

  if ( B_.e_ == 0 )
    B_.e_ = gsl_odeiv_evolve_alloc( State_::STATE_VEC_SIZE );
  else
    gsl_odeiv_evolve_reset( B_.e_ );

  B_.sys_.function = aeif_cond_exp_gridprecise_dynamics;
  B_.sys_.jacobian = NULL;
  B_.sys_.dimension = State_::STATE_VEC_SIZE;
  B_.sys_.params = reinterpret_cast< void* >( this );

  B_.I_stim_ = 0.0;
}

void
nest::aeif_cond_exp_gridprecise::calibrate()
{
  B_.logger_.init(); // ensures initialization in case mm connected after Simulate

  V_.RefractoryCounts_ = Time( Time::ms( P_.t_ref_ ) ).get_steps();
  V_.RefractoryOffset_ = P_.t_ref_ - V_.RefractoryCounts_ * Time::get_resolution().get_ms();
  assert( V_.RefractoryCounts_ >= 0 ); // since t_ref_ >= 0, this can only fail in error
  assert( V_.RefractoryOffset_ >= 0. );
}

/* ----------------------------------------------------------------
 * Update and spike handling functions
 * ---------------------------------------------------------------- */

void
nest::aeif_cond_exp_gridprecise::interpolate_( double& t, double t_old )
{
  // find the exact time when the threshold was crossed
  double dt_crossing = ( P_.V_peak_ - S_.y_old_[ State_::V_M ] ) * ( t - t_old )
    / ( S_.y_[ State_::V_M ] - S_.y_old_[ State_::V_M ] );

  // reset V_m and set the other variables correctly
  S_.y_[ State_::V_M ] = P_.V_reset_;
  for ( int i = 1; i < State_::STATE_VEC_SIZE; ++i )
  {
    S_.y_[ i ] = S_.y_old_[ i ] + ( S_.y_[ i ] - S_.y_old_[ i ] ) / ( t - t_old ) * dt_crossing;
  }
  S_.y_[ State_::W ] += P_.b; // spike-driven adaptation

  t = t_old + dt_crossing;
}

void
nest::aeif_cond_exp_gridprecise::spiking_( Time const& origin, const long_t lag, const double t )
{
  // spike event
  set_spiketime( Time::step( origin.get_steps() + lag + 1 ) );
  SpikeEvent se;
  kernel().event_delivery_manager.send( *this, se, lag );

  // refractoriness
  if ( P_.t_ref_ > 0. )
  {
    S_.r_ = V_.RefractoryCounts_;
    S_.r_offset_ = V_.RefractoryOffset_ - ( B_.step_ - t );
    if ( S_.r_offset_ < 0. )
    {
      if ( S_.r_ > 0 )
      {
        --S_.r_;
        S_.r_offset_ = B_.step_ + S_.r_offset_;
      }
      else
        S_.r_offset_ = t + V_.RefractoryOffset_;
    }
  }
}

void
nest::aeif_cond_exp_gridprecise::update( const Time& origin,
  const nest::long_t from,
  const nest::long_t to )
{
  assert( to >= 0 && ( delay ) from < kernel().connection_builder_manager.get_min_delay() );
  assert( from < to );
  assert( State_::V_M == 0 );

  double t, t_old, t_next_event;

  /* Neurons may have been initialized to superthreshold potentials.
     We need to check for this here and issue spikes at the beginning of
     the interval.
  */
  if ( S_.y_[ State_::V_M ] >= P_.V_peak_ )
  {
    S_.y_[ State_::V_M ] = P_.V_reset_;
    S_.y_[ State_::W ] += P_.b;
    set_spiketime( Time::step( origin.get_steps() + from + 1 ) );
    SpikeEvent se;
    kernel().event_delivery_manager.send( *this, se, from );
  }

  for ( long_t lag = from; lag < to; ++lag )
  {
    t = 0.;

    if ( S_.r_ > 0 )
      --S_.r_;

    // numerical integration with adaptive step size control:
    // ------------------------------------------------------
    // gsl_odeiv_evolve_apply performs only a single numerical
    // integration step, starting from t and bounded by step;
    // the while-loop ensures integration over the whole simulation
    // step (0, step] if more than one integration step is needed due
    // to a small integration step size;
    // note that (t+IntegrationStep > step) leads to integration over
    // (t, step] and afterwards setting t to step, but it does not
    // enforce setting IntegrationStep to step-t

    while ( t < B_.step_ )
    {
      // store the previous values of V_m, w, and t
      std::copy( S_.y_, S_.y_ + sizeof( S_.y_ ) / sizeof( S_.y_[ 0 ] ), S_.y_old_ );
      t_old = t;

      // check for end of refractory period
      if ( P_.t_ref_ > 0. && S_.r_ == 0 && t < S_.r_offset_ )
        t_next_event = S_.r_offset_;
      else
        t_next_event = B_.step_;

      while ( t < t_next_event )
      {
        // propagate the ODE
        const int status = gsl_odeiv_evolve_apply( B_.e_,
          B_.c_,
          B_.s_,
          &B_.sys_,             // system of ODE
          &t,                   // from t
          t_next_event,         // to t <= t_next_event
          &B_.IntegrationStep_, // integration step size
          S_.y_ );              // neuronal state

        // checks
        if ( status != GSL_SUCCESS )
          throw GSLSolverFailure( get_name(), status );
        if ( S_.y_[ State_::V_M ] < -1e3 || S_.y_[ State_::W ] < -1e6 || S_.y_[ State_::W ] > 1e6 )
          throw NumericalInstability( get_name() );
      }

      // check refractoriness
      if ( S_.r_ > 0 || S_.r_offset_ > 0. )
        S_.y_[ State_::V_M ] = P_.V_reset_; // only V_m is frozen
      else if ( S_.y_[ State_::V_M ] >= P_.V_peak_ )
      {
        interpolate_( t, t_old );
        spiking_( origin, lag, t );
      }

      /* reset refractory offset once refractory period is elapsed;
       * this cannot be done beforehand because of the previous check */
      if ( S_.r_ == 0 && std::abs( t - S_.r_offset_ ) < std::numeric_limits< double >::epsilon() )
        S_.r_offset_ = 0.;
    }

    // influence of received spikes on post-synaptic conductances
    S_.y_[ State_::G_EXC ] += B_.spike_exc_.get_value( lag );
    S_.y_[ State_::G_INH ] += B_.spike_inh_.get_value( lag );

    // set new input current
    B_.I_stim_ = B_.currents_.get_value( lag );

    // log state data
    B_.logger_.record_data( origin.get_steps() + lag );
  }
}

void
nest::aeif_cond_exp_gridprecise::handle( SpikeEvent& e )
{
  assert( e.get_delay() > 0 );

  if ( e.get_weight() > 0.0 )
    B_.spike_exc_.add_value(
      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      e.get_weight() * e.get_multiplicity() );
  else
    B_.spike_inh_.add_value(
      e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ),
      -e.get_weight() * e.get_multiplicity() ); // keep conductances positive
}

void
nest::aeif_cond_exp_gridprecise::handle( CurrentEvent& e )
{
  assert( e.get_delay() > 0 );

  const double_t c = e.get_current();
  const double_t w = e.get_weight();

  // add weighted current; HEP 2002-10-04
  B_.currents_.add_value(
    e.get_rel_delivery_steps( kernel().simulation_manager.get_slice_origin() ), w * c );
}

void
nest::aeif_cond_exp_gridprecise::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}

#endif // HAVE_GSL_1_11
