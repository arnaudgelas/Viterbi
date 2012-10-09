#ifndef __viterbi_hpp
#define __viterbi_hpp

#include "viterbi.h"

#include <boost/foreach.hpp>
#include <boost/lexical_cast.hpp>


namespace Viterbi
{
  template< class TState, class TObservation, class TProbability >
  std::ostream&
  operator << (std::ostream& os, HMM< TState, TObservation, TProbability >& h)
  {
    typedef HMM< TState, TObservation, TProbability >             HMMType;
    typedef typename HMMType::StateVectorType                     StateVectorType;
    typedef typename HMMType::ObservationVectorType               ObservationVectorType;
    typedef typename HMMType::StateProbabilityMapType             StateProbabilityMapType;
    typedef typename HMMType::StateStateProbabilityMapType        StateStateProbabilityMapType;
    typedef typename HMMType::ObservationProbabilityMapType       ObservationProbabilityMapType;
    typedef typename HMMType::StateObservationProbabilityMapType  StateObservationProbabilityMapType;

    os << "States:" << std::endl;
    StateVectorType states = h.GetStates();
    BOOST_FOREACH(typename StateVectorType::value_type s, states )
      {
      os << "\tS: " << s << std::endl;
      }

    os << "Observations:" << std::endl;
    ObservationVectorType observations = h.GetObservations();
    BOOST_FOREACH(typename ObservationVectorType::value_type s, observations )
      {
      os << "\tO: " << s << std::endl;
      }

    os << "Start probabilities:" << std::endl;
    StateProbabilityMapType start_probability = h.GetStartProbability();
    BOOST_FOREACH(typename StateProbabilityMapType::value_type m, start_probability)
      {
      os << "\tS: " << m.first
         << " P: " << m.second << std::endl;
      }

    os << "Transition probabilities:" << std::endl;
    StateStateProbabilityMapType transition = h.GetTransitionProbability();

    BOOST_FOREACH(typename StateStateProbabilityMapType::value_type tm, transition)
      {
      BOOST_FOREACH(typename StateProbabilityMapType::value_type m, tm.second)
        {
        os << "\t FS: " << tm.first
           << " TS: " << m.first
           << " P: " << m.second << std::endl;
        }
      }

    os << "Emission probabilities:" << std::endl;
    StateObservationProbabilityMapType emission = h.GetEmissionProbability();
    BOOST_FOREACH(typename StateObservationProbabilityMapType::value_type em, emission)
      {
      BOOST_FOREACH(typename ObservationProbabilityMapType::value_type m, em.second)
        {
        os  << "\tFS: " << em.first
            << " TO: " << m.first
            << " P: " << m.second << std::endl;
        }
      }
    return os;
  }


  //
  // supporting structure
  //
  template< class TState, class TProbability >
  struct Tracking
  {
    typedef TState                    StateType;
    typedef std::vector< StateType >  StateVectorType;

    typedef TProbability              ProbabilityType;

    ProbabilityType prob;
    StateVectorType v_path;
    ProbabilityType v_prob;

    Tracking() : prob(0.0), v_prob(0.0) {}
    Tracking(const ProbabilityType p,
             StateVectorType& v_pth,
             const ProbabilityType v_p):
      prob(p), v_path(v_pth), v_prob(v_p) {}
  };


  // computes total probability for observation
  // most likely viterbi path
  // and probability of such path
  template< class TState, class TObservation, class TProbability >
  void ForwardViberti< TState, TObservation, TProbability >::forward_viterbi()
  {
    m_TotalProbabilityOfObservationSequence = 0;
    m_TotalProbabilityOfVibertiPath         = 0;

    typedef Tracking< TState, TProbability >                TrackingType;
    typedef boost::unordered_map<StateType, TrackingType >  TrackerMapType;
    TrackerMapType T;

    for( typename StateVectorType::iterator sIt = this->m_States.begin();
         sIt != m_States.end();
         ++sIt )
      {
      StateType s = *sIt;

      StateVectorType v_pth( 1 );
      v_pth[0] = s;

      typename StateProbabilityMapType::iterator startIt = m_StartProbability.find( s );

      if( startIt != m_StartProbability.end() )
        {
        ProbabilityType probS = startIt->second;

        T[s] = TrackingType( probS, v_pth, probS );
        }
      else
        {
        std::ostringstream message;
        message << "Error: " << __FILE__ << " " << __LINE__ << std::endl;
        message << "m_StartProbability[" << s <<"] does not exist" << std::endl;

        MyException e( message.str() );
        throw e;
        }
      }

    for( typename ObservationVectorType::iterator obsIt = this->m_Observations.begin();
         obsIt != this->m_Observations.end();
         ++obsIt )
      {
      ObservationType output = *obsIt;
      TrackerMapType U;

      for( typename StateVectorType::iterator sIt = this->m_States.begin();
           sIt != m_States.end();
           ++sIt )
        {
        StateType next_state = *sIt;
        TrackingType next_tracker;

        for( typename StateVectorType::iterator sIt2 = this->m_States.begin();
             sIt2 != m_States.end();
             ++sIt2 )
          {
          StateType source_state = *sIt2;
          TrackingType source_tracker = T[source_state];

          ProbabilityType p;

          typename StateObservationProbabilityMapType::iterator eIt = this->m_EmissionProbability.find( source_state );

          if( eIt != this->m_EmissionProbability.end() )
            {
            typename ObservationProbabilityMapType::iterator oeIt = ( eIt->second ).find( output );

            if( oeIt != ( eIt->second ).end() )
              {
              p = oeIt->second;
              }
            else
              {
              std::ostringstream message;
              message << "Error: " << __FILE__ << " " << __LINE__ << std::endl;
              message << "this->m_EmissionProbability[" << source_state
                      << "][" << output << "] does not exist" << std::endl;

              MyException e( message.str() );
              throw e;
              }
            }
          else
            {
            std::ostringstream message;
            message << "Error: " << __FILE__ << " " << __LINE__ << std::endl;
            message << "this->m_EmissionProbability[" << source_state
                    << "][" << output << "] does not exist" << std::endl;

            MyException e( message.str() );
            throw e;
            }

          typename StateStateProbabilityMapType::iterator tIt = this->m_TransitionProbability.find( source_state );
          if( tIt != this->m_TransitionProbability.end() )
            {
            typename StateProbabilityMapType::iterator stIt = ( tIt->second ).find( next_state );

            if( stIt != ( tIt->second ).end() )
              {
              p *= stIt->second;
              }
            else
              {
              std::ostringstream message;
              message << "Error: " << __FILE__ << " " << __LINE__ << std::endl;
              message << "this->m_TransitionProbability[" << source_state
                      << "][" << next_state << "] does not exist" << std::endl;

              MyException e( message.str() );
              throw e;
              }
            }
          else
            {
            std::ostringstream message;
            message << "Error: " << __FILE__ << " " << __LINE__ << std::endl;
            message << "this->m_TransitionProbability[" << source_state
                    << "][" << next_state << "] does not exist" << std::endl;

            MyException e( message.str() );
            throw e;
            }

          source_tracker.prob *= p;
          source_tracker.v_prob *= p;

          next_tracker.prob += source_tracker.prob;

          if(source_tracker.v_prob > next_tracker.v_prob)
            {
            next_tracker.v_path = source_tracker.v_path;
            next_tracker.v_path.push_back(next_state);
            next_tracker.v_prob = source_tracker.v_prob;
            }
          }
        U[next_state] = next_tracker;
        }

      T = U;
      }

    // apply sum/max to the final states
    TrackingType final_tracker;

    for(typename StateVectorType::iterator sIt = this->m_States.begin();
        sIt != this->m_States.end();
        ++sIt )
      {
      StateType state = *sIt;
      typename TrackerMapType::iterator tIt = T.find( state );

      if( tIt != T.end() )
        {
        TrackingType tracker = tIt->second;

        final_tracker.prob += tracker.prob;

        if( tracker.v_prob > final_tracker.v_prob )
          {
          final_tracker.v_path = tracker.v_path;
          final_tracker.v_prob = tracker.v_prob;
          }
        }
      else
        {
        std::ostringstream message;
        message << "Error: " << __FILE__ << ":" << __LINE__ << std::endl;
        message << "T[" << state << "] does not exist" << std::endl;

        MyException e( message.str() );
        throw e;
        }
      }

    m_TotalProbabilityOfObservationSequence = final_tracker.prob;
    m_TotalProbabilityOfVibertiPath         = final_tracker.v_prob;

    m_FinalPath = StateVectorType( final_tracker.v_path.begin(), final_tracker.v_path.end() );
  }

} // end namespace viterbi
#endif
