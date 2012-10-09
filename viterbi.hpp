#ifndef __viterbi_hpp
#define __viterbi_hpp

#include "viterbi.h"

#include <boost/foreach.hpp>


namespace Viterbi
{

  template< class TState, class TProbability >
  std::ostream&
  operator << (std::ostream& os, HMM< TState, TProbability >& h)
  {
    typedef HMM< TState, TProbability >                     HMMType;
    typedef typename HMMType::StateVectorType               StateVectorType;
    typedef typename HMMType::StateProbabilityMapType       StateProbabilityMapType;
    typedef typename HMMType::StateStateProbabilityMapType  StateStateProbabilityMapType;

    os << "States:" << std::endl;
    StateVectorType states = h.GetStates();
    BOOST_FOREACH(typename StateVectorType::value_type s, states )
      {
      os << "\tS: " << s << std::endl;
      }

    os << "Observations:" << std::endl;
    StateVectorType observations = h.GetObservations();
    BOOST_FOREACH(typename StateVectorType::value_type s, observations )
      {
      os << "\tO: " << s << std::endl;
      }

    os << "Start probabilities:" << std::endl;
    StateProbabilityMapType start_probability = h.GetStartProbability();
    BOOST_FOREACH(typename StateProbabilityMapType::value_type m, start_probability)
      {
      os << "\tS: " << m.first << " P: " << m.second << std::endl;
      }

    os << "Transition probabilities:" << std::endl;
    StateStateProbabilityMapType transition = h.GetTransitionProbability();

    BOOST_FOREACH(typename StateStateProbabilityMapType::value_type tm, transition)
      {
      BOOST_FOREACH(typename StateProbabilityMapType::value_type m, tm.second)
        {
        os << "\t FS: " << tm.first << " TS: " << m.first
           << " P: " << m.second << std::endl;
        }
      }

    os << "Emission probabilities:" << std::endl;
    StateStateProbabilityMapType emission = h.GetEmissionProbability();
    BOOST_FOREACH(typename StateStateProbabilityMapType::value_type em, emission)
      {
      BOOST_FOREACH(typename StateProbabilityMapType::value_type m,
          em.second)
        {
        os  << "\tFS: " << em.first << " TO: " << m.first
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
  template< class TState, class TProbability >
  void ForwardViberti< TState, TProbability >::forward_viterbi()
  {
    m_TotalProbabilityOfObservationSequence = 0;
    m_TotalProbabilityOfVibertiPath         = 0;

    typedef Tracking< TState, TProbability >                TrackingType;
    typedef boost::unordered_map<StateType, TrackingType >  TrackerMapType;
    TrackerMapType T;

    BOOST_FOREACH(typename StateVectorType::value_type s, this->m_States)
      {
      StateVectorType v_pth( 1 );
      v_pth[0] = s;

      ProbabilityType probS = m_StartProbability[s];

      T[s] = TrackingType( probS, v_pth, probS );
      }

    BOOST_FOREACH(typename StateVectorType::value_type output, this->m_Observations)
      {
      TrackerMapType U;

      BOOST_FOREACH(typename StateVectorType::value_type next_state, this->m_States)
        {
        TrackingType next_tracker;

        BOOST_FOREACH(typename StateVectorType::value_type source_state, this->m_States)
          {
          TrackingType source_tracker = T[source_state];

          ProbabilityType p = this->m_EmissionProbability[source_state][output];
          p *= this->m_TransitionProbability[source_state][next_state];

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

    BOOST_FOREACH(typename StateVectorType::value_type state,  this->m_States)
      {
      TrackingType tracker = T[state];

      final_tracker.prob += tracker.prob;

      if(tracker.v_prob > final_tracker.v_prob)
        {
        final_tracker.v_path = tracker.v_path;
        final_tracker.v_prob = tracker.v_prob;
        }
      }

    m_TotalProbabilityOfObservationSequence = final_tracker.prob;
    m_TotalProbabilityOfVibertiPath         = final_tracker.v_prob;

    m_FinalPath = StateVectorType( final_tracker.v_path.begin(), final_tracker.v_path.end() );
  }

} // end namespace viterbi
#endif
