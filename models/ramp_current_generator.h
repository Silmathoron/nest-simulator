/*
 *  ramp_current_generator.h
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


#ifndef RAMP_CURRENT_GENERATOR_H
#define RAMP_CURRENT_GENERATOR_H

// C++ includes:
#include <vector>

// Includes from nestkernel:
#include "connection.h"
#include "device_node.h"
#include "event.h"
#include "nest_types.h"
#include "ring_buffer.h"
#include "stimulating_device.h"
#include "universal_data_logger.h"

namespace nest
{

/** @BeginDocumentation
@ingroup Devices
@ingroup generator

Name: ramp_current_generator - provides a ramp input current

Description:

The ramp_current_generator provides a ramp current input, i.e. a current which
varies linearly over time, to the connected node(s).
The unit of the current is pA.

Parameters:

The following parameters can be set in the status dictionary:

\verbatim embed:rst
==================== ===============  ========================================
 amplitude_times     list of ms       Times at which current changes
 amplitude_values    list of pA       Amplitudes of the current at these times
==================== ===============  ========================================
\endverbatim

Notes:

Times of amplitude changes must be strictly increasing.
An initial value of 0 pA at time 0 ms is always considered and, if used,
does not need to be specified.

Examples:

The current can be configured in the following way:

    import nest

    resol = 0.1
    nest.SetKernelStatus({"resolution": resol})

    nest.Create("ramp_current_generator", params={
        "amplitude_times": [100., 100. + resol, 110., 140.],
        "amplitude_values": [10., 0., 0., -5.],
        "stop": 140. + resol
    })

The amplitude of the current will increases from 0 to 10 pA in 100 ms, then
go down to zero and remain there for 10 ms, then go from 0 to -5 pA in
30 ms before going back to zero until the end of the simulation.

Sends: CurrentEvent

Author: Tanguy Fardet

SeeAlso: ac_generator, dc_generator, step_current_generator, Device,
StimulatingDevice
*/
class ramp_current_generator : public DeviceNode
{

public:
  ramp_current_generator();
  ramp_current_generator( const ramp_current_generator& );

  bool
  has_proxies() const
  {
    return false;
  }

  port send_test_event( Node&, rport, synindex, bool );

  using Node::handle;
  using Node::handles_test_event;

  void handle( DataLoggingRequest& );

  port handles_test_event( DataLoggingRequest&, rport );

  void get_status( DictionaryDatum& ) const;
  void set_status( const DictionaryDatum& );

  //! Allow multimeter to connect to local instances
  bool
  local_receiver() const
  {
    return true;
  }

private:
  void init_state_( const Node& );
  void init_buffers_();
  void calibrate();

  void update( Time const&, const long, const long );

  struct Buffers_;

  /**
   * Store independent parameters of the model.
   */
  struct Parameters_
  {
    //! Times of amplitude changes
    std::vector< Time > amp_time_stamps_;

    //! Amplitude values activated at given times
    std::vector< double > amp_values_;

    Parameters_(); //!< Sets default parameter values
    Parameters_( const Parameters_&, Buffers_& );
    Parameters_( const Parameters_& );
    Parameters_& operator=( const Parameters_& p );

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
    //! Set values from dictionary
    void set( const DictionaryDatum&, Buffers_& );

    /**
     * Return time as Time object if valid, otherwise throw BadProperty
     *
     * @param amplitude time, ms
     * @param previous time stamp
     */
    Time validate_time_( double, const Time& );
  };

  // ------------------------------------------------------------

  struct State_
  {
    double I_; //!< Instantaneous current value; used for recording current

    State_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
  };

  // ------------------------------------------------------------

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< ramp_current_generator >;
  friend class UniversalDataLogger< ramp_current_generator >;

  // ------------------------------------------------------------

  struct Buffers_
  {
    size_t idx_;   //!< index of current amplitude
    long stp_;     //!< step associated to initial amplitude
    double amp_;   //!< current initial amplitude
    double slope_; //!< current slope value (pA/step)

    Buffers_( ramp_current_generator& );
    Buffers_( const Buffers_&, ramp_current_generator& );
    UniversalDataLogger< ramp_current_generator > logger_;
  };

  // ------------------------------------------------------------

  double
  get_I_() const
  {
    return S_.I_;
  }

  // ------------------------------------------------------------

  StimulatingDevice< CurrentEvent > device_;
  static RecordablesMap< ramp_current_generator > recordablesMap_;
  Parameters_ P_;
  State_ S_;
  Buffers_ B_;
};

inline port
ramp_current_generator::send_test_event( Node& target, rport receptor_type, synindex syn_id, bool )
{
  device_.enforce_single_syn_type( syn_id );

  CurrentEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
ramp_current_generator::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
ramp_current_generator::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  device_.get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
ramp_current_generator::set_status( const DictionaryDatum& d )
{
  Parameters_ ptmp = P_; // temporary copy in case of errors
  ptmp.set( d, B_ );     // throws if BadProperty

  // We now know that ptmp is consistent. We do not write it back
  // to P_ before we are also sure that the properties to be set
  // in the parent class are internally consistent.
  device_.set_status( d );

  // if we get here, temporaries contain consistent set of properties
  P_ = ptmp;
}

} // namespace

#endif /* #ifndef ramp_current_generator_H */
