/*
 *  ramp_current_generator.cpp
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

#include "ramp_current_generator.h"

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
RecordablesMap< ramp_current_generator > ramp_current_generator::recordablesMap_;

template <>
void
RecordablesMap< ramp_current_generator >::create()
{
  insert_( Name( names::I ), &ramp_current_generator::get_I_ );
}
}

/* ----------------------------------------------------------------
 * Default constructors defining default parameter
 * ---------------------------------------------------------------- */

nest::ramp_current_generator::Parameters_::Parameters_()
  : amp_time_stamps_()
  , amp_values_() // pA
{
}

nest::ramp_current_generator::Parameters_::Parameters_( const Parameters_& p )
  : amp_time_stamps_( p.amp_time_stamps_ )
  , amp_values_( p.amp_values_ )
{
}

nest::ramp_current_generator::Parameters_& nest::ramp_current_generator::Parameters_::operator=( const Parameters_& p )
{
  if ( this == &p )
  {
    return *this;
  }

  amp_time_stamps_ = p.amp_time_stamps_;
  amp_values_ = p.amp_values_;

  return *this;
}

nest::ramp_current_generator::State_::State_()
  : I_( 0.0 ) // pA
{
}

nest::ramp_current_generator::Buffers_::Buffers_( ramp_current_generator& n )
  : idx_( 0 )
  , stp_( 0 )
  , amp_( 0.0 )
  , slope_( 0.0 ) // pA/step
  , logger_( n )
{
}

nest::ramp_current_generator::Buffers_::Buffers_( const Buffers_&, ramp_current_generator& n )
  : idx_( 0 )
  , stp_( 0 )
  , amp_( 0.0 )
  , slope_( 0.0 ) // pA/step
  , logger_( n )
{
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */

void
nest::ramp_current_generator::Parameters_::get( DictionaryDatum& d ) const
{
  std::vector< double >* times_ms = new std::vector< double >();
  times_ms->reserve( amp_time_stamps_.size() );
  for ( std::vector< Time >::const_iterator it = amp_time_stamps_.begin(); it != amp_time_stamps_.end(); ++it )
  {
    times_ms->push_back( it->get_ms() );
  }
  ( *d )[ names::amplitude_times ] = DoubleVectorDatum( times_ms );
  ( *d )[ names::amplitude_values ] = DoubleVectorDatum( new std::vector< double >( amp_values_ ) );
}

nest::Time
nest::ramp_current_generator::Parameters_::validate_time_( double t, const Time& t_previous )
{
  if ( t <= 0.0 )
  {
    throw BadProperty(
      "Amplitude can only be changed at strictly "
      "positive times (t > 0)." );
  }

  // Force the amplitude change time to the grid
  // First, convert the time to tics, may not be on grid
  Time t_amp = Time::ms( t );
  if ( not t_amp.is_grid_time() )
  {
      std::stringstream msg;
      msg << "ramp_current_generator: Time point " << t << " is not representable in current resolution.";
      throw BadProperty( msg.str() );
  }

  assert( t_amp.is_grid_time() );

  // t_amp is now the correct time stamp given the chosen options
  if ( t_amp <= t_previous )
  {
    throw BadProperty(
      "ramp_current_generator: amplitude "
      "times must be at strictly increasing "
      "time steps." );
  }

  // when we get here, we know that the spike time is valid
  return t_amp;
}

void
nest::ramp_current_generator::Parameters_::set( const DictionaryDatum& d, Buffers_& b )
{
  std::vector< double > new_times;
  const bool times_changed = updateValue< std::vector< double > >( d, names::amplitude_times, new_times );
  const bool values_changed = updateValue< std::vector< double > >( d, names::amplitude_values, amp_values_ );

  if ( times_changed xor values_changed )
  {
    throw BadProperty( "Amplitude times and values must be reset together." );
  }

  const size_t times_size = times_changed ? new_times.size() : amp_time_stamps_.size();

  if ( times_size != amp_values_.size() )
  {
    throw BadProperty( "Amplitude times and values have to be the same size." );
  }

  if ( times_changed )
  {
    std::vector< Time > new_stamps;
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
    amp_time_stamps_.swap( new_stamps );
  }

  if ( times_changed or values_changed )
  {
    b.idx_ = 0; // reset if we got new data
  }
}


/* ----------------------------------------------------------------
 * Default and copy constructor for node
 * ---------------------------------------------------------------- */

nest::ramp_current_generator::ramp_current_generator()
  : DeviceNode()
  , device_()
  , P_()
  , S_()
  , B_( *this )
{
  recordablesMap_.create();
}

nest::ramp_current_generator::ramp_current_generator( const ramp_current_generator& n )
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
nest::ramp_current_generator::init_state_( const Node& proto )
{
  const ramp_current_generator& pr = downcast< ramp_current_generator >( proto );

  device_.init_state( pr.device_ );
}

void
nest::ramp_current_generator::init_buffers_()
{
  device_.init_buffers();
  B_.logger_.reset();

  B_.idx_   = 0;
  B_.stp_   = 0.;
  B_.amp_   = 0.;
  B_.slope_ = 0.;
}

void
nest::ramp_current_generator::calibrate()
{
  B_.logger_.init();

  device_.calibrate();
}


/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */

void
nest::ramp_current_generator::update( Time const& origin, const long from, const long to )
{
  assert( to >= 0 && ( delay ) from < kernel().connection_manager.get_min_delay() );
  assert( from < to );

  assert( P_.amp_time_stamps_.size() == P_.amp_values_.size() );

  const long t0 = origin.get_steps();

  // Skip any times in the past. Since we must send events proactively,
  // idx_ must point to times in the future.
  const long first = t0 + from;
  while ( B_.idx_ < P_.amp_time_stamps_.size() && P_.amp_time_stamps_[ B_.idx_ ].get_steps() <= first )
  {
    ++B_.idx_;
  }

  for ( long offs = from; offs < to; ++offs )
  {
    const long curr_time = t0 + offs;

    S_.I_ = 0.0;

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
                    / (P_.amp_time_stamps_[ B_.idx_ + 1] - P_.amp_time_stamps_[ B_.idx_ ]);
      }
      else if ( B_.idx_ > 1)
      {
        B_.amp_   = P_.amp_values_[ B_.idx_ - 1 ];
        B_.stp_   = P_.amp_time_stamps_[ B_.idx_ - 1 ].get_steps();
        B_.slope_ = (P_.amp_values_[ B_.idx_ ] - P_.amp_values_[ B_.idx_ - 1 ])
                    / (P_.amp_time_stamps_[ B_.idx_ ] - P_.amp_time_stamps_[ B_.idx_ - 1 ]);
      }
      else
      {
        B_.amp_   = 0.;
        B_.stp_   = 0;
        B_.slope_ = P_.amp_values_[ B_.idx_ ] / P_.amp_time_stamps_[ B_.idx_ ];
      }
      B_.idx_++;
    }

    // but send only if active
    if ( device_.is_active( Time::step( curr_time ) ) )
    {
      double instantaneous_current = B_.amp_ + B_.slope_*(curr_time - B_.stp_);
      CurrentEvent ce;
      ce.set_current( instantaneous_current );
      S_.I_ = instantaneous_current;
      kernel().event_delivery_manager.send( *this, ce, offs );
    }
    B_.logger_.record_data( origin.get_steps() + offs );
  }
}

void
nest::ramp_current_generator::handle( DataLoggingRequest& e )
{
  B_.logger_.handle( e );
}