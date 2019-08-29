/*
 *  sinusoidal_current_generator.cpp
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

#include "sinusoidal_current_generator.h"

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "kernel_manager.h"
#include "universal_data_logger_impl.h"

// Includes from sli:
#include "booldatum.h"
#include "dict.h"
#include "dictutils.h"
#include "doubledatum.h"
#include "integerdatum.h"

namespace nest
{
RecordablesMap< sinusoidal_current_generator > sinusoidal_current_generator::recordablesMap_;

template <>
void
RecordablesMap< sinusoidal_current_generator >::create()
{
  insert_( Name( names::I ), &sinusoidal_current_generator::get_I_ );
  insert_( Name( names::frequency ), &sinusoidal_current_generator::get_f_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

nest::sinusoidal_current_generator::Parameters_::Parameters_()
  : amp_time_stamps_()
  , freq_time_stamps_()
  , phase_time_stamps_()
  , amp_values_() // pA
  , freq_values_() // Hz
  , phase_values_() // real
  , log_sweep_( false )
{
}

nest::sinusoidal_current_generator::Parameters_::Parameters_( const Parameters_& p )
  : amp_time_stamps_( p.amp_time_stamps_ )
  , freq_time_stamps_( p.freq_time_stamps_ )
  , phase_time_stamps_( p.phase_time_stamps_ )
  , amp_values_( p.amp_values_ )
  , freq_values_( p.freq_values_ )
  , phase_values_( p.phase_values_ )
  , log_sweep_( p.log_sweep_ )
{
}

nest::sinusoidal_current_generator::Parameters_& nest::sinusoidal_current_generator::Parameters_::operator=( const Parameters_& p )
{
  if ( this == &p )
  {
    return *this;
  }

  amp_time_stamps_ = p.amp_time_stamps_;
  amp_values_ = p.amp_values_;

  freq_time_stamps_ = p.freq_time_stamps_;
  freq_values_ = p.freq_values_;

  phase_time_stamps_ = p.phase_time_stamps_;
  phase_values_ = p.phase_values_;

  log_sweep_ = p.log_sweep_;

  return *this;
}

nest::sinusoidal_current_generator::State_::State_()
  : I_( 0.0 ) // pA
  , f_( 0.0 ) // Hz
{
}

nest::sinusoidal_current_generator::Buffers_::Buffers_( sinusoidal_current_generator& n )
  : idx_( 0 )
  , stp_( 0 )
  , amp_( 0.0 )   // pA
  , slope_( 0.0 ) // pA/step
  , fidx_( 0 )
  , fstp_( 0 )
  , freq_( 0. )   // Hz
  , fslope_( 0. ) // Hz/step
  , pidx_( 0 )
  , phase_( 0. )  // real
  , logger_( n )
{
}

nest::sinusoidal_current_generator::Buffers_::Buffers_( const Buffers_&, sinusoidal_current_generator& n )
  : idx_( 0 )
  , stp_( 0 )
  , amp_( 0.0 )
  , slope_( 0.0 ) // pA/step
  , fidx_( 0 )
  , fstp_( 0 )
  , freq_( 0. )   // Hz
  , fslope_( 0. ) // Hz/step
  , pidx_( 0 )
  , phase_( 0. )  // real
  , logger_( n )
{
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::sinusoidal_current_generator::Parameters_::get( DictionaryDatum& d ) const
{
  std::vector< double >* atimes_ms = new std::vector< double >();
  atimes_ms->reserve( amp_time_stamps_.size() );

  std::vector< double >* ftimes_ms = new std::vector< double >();
  ftimes_ms->reserve( freq_time_stamps_.size() );

  std::vector< double >* ptimes_ms = new std::vector< double >();
  ptimes_ms->reserve( phase_time_stamps_.size() );

  for ( auto it : amp_time_stamps_ )
  {
    atimes_ms->push_back( it.get_ms() );
  }

  for ( auto it : freq_time_stamps_ )
  {
    ftimes_ms->push_back( it.get_ms() );
  }

  for ( auto it : phase_time_stamps_ )
  {
    ptimes_ms->push_back( it.get_ms() );
  }

  ( *d )[ names::amplitude_times ] = DoubleVectorDatum( atimes_ms );
  ( *d )[ names::amplitude_values ] = DoubleVectorDatum( new std::vector< double >( amp_values_ ) );
  ( *d )[ names::frequency_times ] = DoubleVectorDatum( ftimes_ms );
  ( *d )[ names::frequency_values ] = DoubleVectorDatum( new std::vector< double >( freq_values_ ) );
  ( *d )[ names::phase_times ] = DoubleVectorDatum( ptimes_ms );
  ( *d )[ names::phase_values ] = DoubleVectorDatum( new std::vector< double >( phase_values_ ) );

  ( *d )[ "log_sweep" ] = BoolDatum( log_sweep_ );
}

nest::Time
nest::sinusoidal_current_generator::Parameters_::validate_time_( double t, const Time& t_previous )
{
  if ( t <= 0.0 )
  {
    throw BadProperty(
      "Amplitude can only be changed at strictly "
      "positive times (t > 0)." );
  }

  // Force the amplitude change time to the grid
  // First, convert the time to tics, may not be on grid
  Time t_val = Time::ms( t );
  if ( not t_val.is_grid_time() )
  {
      std::stringstream msg;
      msg << "sinusoidal_current_generator: Time point " << t << " is not representable in current resolution.";
      throw BadProperty( msg.str() );
  }

  assert( t_val.is_grid_time() );

  // t_amp is now the correct time stamp given the chosen options
  if ( t_val <= t_previous )
  {
    throw BadProperty(
      "sinusoidal_current_generator: amplitude "
      "times must be at strictly increasing "
      "time steps." );
  }

  // when we get here, we know that the spike time is valid
  return t_val;
}

void nest::sinusoidal_current_generator::Parameters_::check_and_update_params(
  bool times_changed, bool values_changed, const std::vector< double >& new_times,
  std::vector< Time >& timestamps_to_update, const std::vector< double >& new_values,
  const std::string& name, size_t& idx )
{
  if ( times_changed xor values_changed )
  {
    throw nest::BadProperty( name + " times and values must be reset together." );
  }

  const size_t times_size = times_changed ? new_times.size() : timestamps_to_update.size();

  if ( times_size != new_values.size() )
  {
    throw nest::BadProperty( name + " times and values have to be the same size." );
  }

  if ( times_changed )
  {
    std::vector< nest::Time > new_stamps;
    new_stamps.reserve( times_size );

    if ( not new_times.empty() )
    {
      // insert first change, we are sure we have one
      new_stamps.push_back( validate_time_( new_times[ 0 ], Time( Time::ms( 0 ) ) ) );

      // insert all others
      for ( size_t idx = 1; idx < times_size; ++idx )
      {
        new_stamps.push_back( validate_time_( new_times[ idx ], new_stamps[ idx - 1 ] ) );
      }
    }

    // if we get here, all times have been successfully converted
    timestamps_to_update.swap( new_stamps );
  }

  if ( times_changed or values_changed )
  {
    idx = 0; // reset if we got new data
  }
}

void
nest::sinusoidal_current_generator::Parameters_::set( const DictionaryDatum& d, Buffers_& b )
{
  std::vector< double > new_atimes, new_ftimes, new_ptimes;
  const bool atimes_changed = updateValue< std::vector< double > >( d, names::amplitude_times, new_atimes );
  const bool avalues_changed = updateValue< std::vector< double > >( d, names::amplitude_values, amp_values_ );
  const bool ftimes_changed = updateValue< std::vector< double > >( d, names::frequency_times, new_ftimes );
  const bool fvalues_changed = updateValue< std::vector< double > >( d, names::frequency_values, freq_values_ );
  const bool ptimes_changed = updateValue< std::vector< double > >( d, names::phase_times, new_ptimes );
  const bool pvalues_changed = updateValue< std::vector< double > >( d, names::phase_values, phase_values_ );

  check_and_update_params(atimes_changed, avalues_changed, new_atimes, amp_time_stamps_, amp_values_, "Amplitude",
                          b.idx_);

  check_and_update_params(ftimes_changed, fvalues_changed, new_ftimes, freq_time_stamps_, freq_values_, "Frequency",
                          b.fidx_);

  check_and_update_params(ptimes_changed, pvalues_changed, new_ptimes, phase_time_stamps_, phase_values_, "Phase",
                          b.pidx_);

  updateValue< bool >( d, "log_sweep", log_sweep_ );
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::sinusoidal_current_generator::sinusoidal_current_generator()
  : DeviceNode()
  , device_()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::sinusoidal_current_generator::sinusoidal_current_generator( const sinusoidal_current_generator& n )
  : DeviceNode( n )
  , device_( n.device_ )
  , P_( n.P_ )
  , S_( n.S_ )
  , B_( n.B_, *this )
{
}


/* ----------------------------------------------------------------
 * Node initialization functions
 * ---------------------------------------------------------------- */

void
nest::sinusoidal_current_generator::init_state_( const Node& proto )
{
  const sinusoidal_current_generator& pr = downcast< sinusoidal_current_generator >( proto );

  device_.init_state( pr.device_ );
}

void
nest::sinusoidal_current_generator::init_buffers_()
{
  device_.init_buffers();
  B_.logger_.reset();

  B_.idx_   = 0;
  B_.stp_   = 0.;
  B_.amp_   = 0.;
  B_.slope_ = 0.;

  B_.fidx_   = 0;
  B_.fstp_   = 0.;
  B_.freq_   = 0.;
  B_.fslope_ = 0.;

  B_.pidx_   = 0;
  B_.phase_   = 0.;
}

void
nest::sinusoidal_current_generator::calibrate()
{
  B_.logger_.init();

  device_.calibrate();
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

void
nest::sinusoidal_current_generator::update( Time const& origin, const long from, const long to )
{
  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  assert( P_.amp_time_stamps_.size() == P_.amp_values_.size() );
  assert( P_.freq_time_stamps_.size() == P_.freq_values_.size() );
  assert( P_.phase_time_stamps_.size() == P_.phase_values_.size() );

  const long t0 = origin.get_steps();

  // ------------------ //
  // Select all indices //
  // ------------------ //

  // Skip any times in the past. Since we must send events proactively,
  // idx_ must point to times in the future.
  const long first = t0 + from;
  while ( B_.idx_ < P_.amp_time_stamps_.size() && P_.amp_time_stamps_[ B_.idx_ ].get_steps() <= first )
  {
    ++B_.idx_;
  }

  while ( B_.idx_ < P_.amp_time_stamps_.size() && P_.amp_time_stamps_[ B_.idx_ ].get_steps() <= first )
  {
    ++B_.fidx_;
  }

  while ( B_.idx_ < P_.amp_time_stamps_.size() && P_.amp_time_stamps_[ B_.idx_ ].get_steps() <= first )
  {
    ++B_.pidx_;
  }

  // --------------- //
  // Set the current //
  // --------------- //

  for ( long offs = from; offs < to; ++offs )
  {
    const long curr_time = t0 + offs;

    S_.I_ = 0.0;

    // --------------- //
    // Amplitudes step //
    // --------------- //

    // Keep the amplitude up-to-date at all times.
    // We need to change the amplitude one step ahead of time, see comment
    // on class SimulatingDevice.
    if ( B_.idx_ < P_.amp_time_stamps_.size() && curr_time + 1 == P_.amp_time_stamps_[ B_.idx_ ].get_steps() )
    {
      if (B_.idx_ < P_.amp_time_stamps_.size() - 1)
      {
        B_.amp_   = P_.amp_values_[ B_.idx_ ];
        B_.stp_   = P_.amp_time_stamps_[ B_.idx_ ].get_steps();
        B_.slope_ = (P_.amp_values_[ B_.idx_ + 1] - P_.amp_values_[ B_.idx_ ])
                    / (P_.amp_time_stamps_[ B_.idx_ + 1].get_steps() - P_.amp_time_stamps_[ B_.idx_ ].get_steps());
      }
      else if ( B_.idx_ > 1)
      {
        B_.amp_   = P_.amp_values_[ B_.idx_ - 1 ];
        B_.stp_   = P_.amp_time_stamps_[ B_.idx_ - 1 ].get_steps();
        B_.slope_ = (P_.amp_values_[ B_.idx_ ] - P_.amp_values_[ B_.idx_ - 1 ])
                    / (P_.amp_time_stamps_[ B_.idx_ ].get_steps() - P_.amp_time_stamps_[ B_.idx_ - 1 ].get_steps());
      }
      else
      {
        B_.amp_   = 0.;
        B_.stp_   = 0;
        B_.slope_ = P_.amp_values_[ B_.idx_ ] / P_.amp_time_stamps_[ B_.idx_ ].get_steps();
      }
      B_.idx_++;
    }

    // -------------- //
    // Frequency step //
    // -------------- //

    // Keep the frequency up-to-date at all times.
    // We need to change the frequency one step ahead of time, see comment
    // on class SimulatingDevice.
    if ( B_.fidx_ < P_.freq_time_stamps_.size() && curr_time + 1 == P_.freq_time_stamps_[ B_.fidx_ ].get_steps() )
    {
      if (B_.fidx_ < P_.freq_time_stamps_.size() - 1)
      {
        B_.freq_   = P_.freq_values_[ B_.fidx_ ];
        B_.fstp_   = P_.freq_time_stamps_[ B_.fidx_ ].get_steps();

        if (P_.log_sweep_)
        {
          B_.fslope_ = (P_.freq_values_[ B_.fidx_ + 1] - B_.freq_)
                      / (std::exp(P_.freq_time_stamps_[ B_.fidx_ + 1].get_steps()) - std::exp(B_.fstp_));
        }
        else
        {
          B_.fslope_ = (P_.freq_values_[ B_.fidx_ + 1] - B_.freq_)
                       / (P_.freq_time_stamps_[ B_.fidx_ + 1].get_steps() - B_.fstp_);
        }
      }
      else if ( B_.fidx_ > 1)
      {
        B_.freq_   = P_.freq_values_[ B_.fidx_ - 1 ];
        B_.fstp_   = P_.freq_time_stamps_[ B_.fidx_ - 1 ].get_steps();

        if (P_.log_sweep_)
        {
          B_.fslope_ = (P_.freq_values_[ B_.fidx_ ] - B_.freq_)
                      / (std::exp(P_.freq_time_stamps_[ B_.fidx_ ].get_steps()) - std::exp(B_.fstp_));
        }
        else
        {
          B_.fslope_ = (P_.freq_values_[ B_.fidx_ ] - B_.freq_)
                      / (P_.freq_time_stamps_[ B_.fidx_ ].get_steps() - B_.fstp_);
        }
      }
      else
      {
        B_.freq_   = 0.;
        B_.fstp_   = 0;

        if (P_.log_sweep_)
        {
          B_.fslope_ = P_.freq_values_[ B_.fidx_ ] / std::exp(P_.freq_time_stamps_[ B_.fidx_ ].get_steps());
        }
        else
        {
          B_.fslope_ = P_.freq_values_[ B_.fidx_ ] / P_.freq_time_stamps_[ B_.fidx_ ].get_steps();
        }
      }
      B_.fidx_++;
    }

    // ---------- //
    // Phase step //
    // ---------- //

    // Keep the frequency up-to-date at all times.
    // We need to change the frequency one step ahead of time, see comment
    // on class SimulatingDevice.
    if ( B_.pidx_ < P_.phase_time_stamps_.size() && curr_time + 1 == P_.phase_time_stamps_[ B_.pidx_ ].get_steps() )
    {
      B_.phase_   = P_.phase_values_[ B_.pidx_ ];
      B_.pidx_++;
    }

    // but send only if active
    if ( device_.is_active( Time::step( curr_time ) ) )
    {
      double amplitude = B_.amp_ + B_.slope_*(curr_time - B_.stp_ + 1);

      double frequency = P_.log_sweep_
        ? B_.freq_ + B_.fslope_*(std::exp(curr_time + 1) - std::exp(B_.fstp_))
        : B_.freq_ + B_.fslope_*(curr_time - B_.fstp_ + 1);

      double time_s    = 0.001*curr_time*Time::get_resolution().get_ms();
      double current   = amplitude*std::sin(frequency*time_s + B_.phase_);

      CurrentEvent ce;
      ce.set_current( current );
      S_.I_ = current;
      kernel().event_delivery_manager.send( *this, ce, offs );
    }
    B_.logger_.record_data( origin.get_steps() + offs );
  }
}

void
nest::sinusoidal_current_generator::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}
