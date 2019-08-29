/*
 *  sinusoidal_current_generator.h
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


#ifndef SINUSOIDAL_CURRENT_GENERATOR_H
#define SINUSOIDAL_CURRENT_GENERATOR_H

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

Name: sinusoidal_current_generator - provides a sinusoidal input current

Description:

The sinusoidal_current_generator provides a sinusoidal current input, i.e. a current which
varies as a sine of time, to the connected node(s).
The frequency of the current can vary linearly or logarithmically with time.
The unit of the current is pA, that of the frequency is Hz.
Both current and frequency can be recorded by a multimeter, under the respective names "I" and "frequency".

Parameters:

The following parameters can be set in the status dictionary:

\verbatim embed:rst
==================== ===============  ========================================
 amplitude_times     list of ms       Times at which current amplitude changes
 amplitude_values    list of pA       Amplitudes of the sinusoidal current
 frequency_times     list of ms       Times at which current frequency changes
 frequency_values    list of Hz       Frequency of the sinusoidal current
 phase_times         list of ms       Times at which current phase changes
 phase_values        list of real     Phase of the sinusoidal current
 log_sweep           bool             Change frequency logarithmically
==================== ===============  ========================================
\endverbatim

Notes:

Times of amplitude, frequency, and phase changes must be strictly increasing.
An initial value of 0 pA, 0 Hz, and 0 phase at time 0 ms is always considered
and, does not need to be specified.
By default, ``log_sweep`` is set to ``False``.

Examples:

The current can be altered in the following way:

    /sinusoidal_current_generator Create /sc Set
    sc << /amplitude_times [0.2 0.5] /amplitude_values [2.0 4.0] >> SetStatus

    The amplitude of the DC will be 0.0 pA in the time interval [0, 0.2),
    2.0 pA in the interval [0.2, 0.5) and 4.0 from then on.

Sends: CurrentEvent

Author: Tanguy Fardet

SeeAlso: ac_generator, dc_generator, step_current_generator, Device,
StimulatingDevice
*/
class sinusoidal_current_generator : public DeviceNode
{

public:
  sinusoidal_current_generator();
  sinusoidal_current_generator( const sinusoidal_current_generator& );

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
    //! Times of amplitude, frequency, and phase changes
    std::vector< Time > amp_time_stamps_;
    std::vector< Time > freq_time_stamps_;
    std::vector< Time > phase_time_stamps_;

    //! Amplitude, frequency, and phase values activated at given times
    std::vector< double > amp_values_;
    std::vector< double > freq_values_;
    std::vector< double > phase_values_;

    bool log_sweep_;  //!< sweep mode ("lin" if false, "log" if true)

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

    /**
     * Check and update the parameters of the generator
     */
    void check_and_update_params(
      bool times_changed, bool values_changed, const std::vector< double > &new_times,
      std::vector< Time > &timestamps_to_update, const std::vector< double > &new_values,
      const std::string &name, size_t& idx );
  };

  // ------------------------------------------------------------

  struct State_
  {
    double I_; //!< Instantaneous current value; used for recording current
    double f_; //!< Instantaneous current value; used for recording current

    State_(); //!< Sets default parameter values

    void get( DictionaryDatum& ) const; //!< Store current values in dictionary
  };

  // ------------------------------------------------------------

  // The next two classes need to be friends to access the State_ class/member
  friend class RecordablesMap< sinusoidal_current_generator >;
  friend class UniversalDataLogger< sinusoidal_current_generator >;

  // ------------------------------------------------------------

  struct Buffers_
  {
    size_t idx_;    //!< index of current amplitude
    long stp_;      //!< step associated to initial amplitude
    double amp_;    //!< current initial amplitude (pA)
    double slope_;  //!< current slope value (pA/step)
    size_t fidx_;   //!< index of current initial frequency
    long fstp_;     //!< step associated to current initial frequency
    double freq_;   //!< current initial frequency (Hz)
    double fslope_; //!< frequency slope value (Hz/step)
    size_t pidx_;   //!< index of current phase
    double phase_;  //!< current phase (real)

    Buffers_( sinusoidal_current_generator& );
    Buffers_( const Buffers_&, sinusoidal_current_generator& );
    UniversalDataLogger< sinusoidal_current_generator > logger_;
  };

  // ------------------------------------------------------------

  double
  get_I_() const
  {
    return S_.I_;
  }

  double
  get_f_() const
  {
    return S_.f_;
  }

  // ------------------------------------------------------------

  StimulatingDevice< CurrentEvent > device_;
  static RecordablesMap< sinusoidal_current_generator > recordablesMap_;
  Parameters_ P_;
  State_ S_;
  Buffers_ B_;
};

inline port
sinusoidal_current_generator::send_test_event( Node& target, rport receptor_type, synindex syn_id, bool )
{
  device_.enforce_single_syn_type( syn_id );

  CurrentEvent e;
  e.set_sender( *this );

  return target.handles_test_event( e, receptor_type );
}

inline port
sinusoidal_current_generator::handles_test_event( DataLoggingRequest& dlr, rport receptor_type )
{
  if ( receptor_type != 0 )
  {
    throw UnknownReceptorType( receptor_type, get_name() );
  }
  return B_.logger_.connect_logging_device( dlr, recordablesMap_ );
}

inline void
sinusoidal_current_generator::get_status( DictionaryDatum& d ) const
{
  P_.get( d );
  device_.get_status( d );

  ( *d )[ names::recordables ] = recordablesMap_.get_list();
}

inline void
sinusoidal_current_generator::set_status( const DictionaryDatum& d )
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

#endif /* #ifndef sinusoidal_current_generator_H */
