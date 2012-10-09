#ifndef __viberti_h
#define __viberti_h

/*
  Implementation of http://en.wikipedia.org/wiki/Viterbi_algorithm
   It simplifies http://bozskyfilip.blogspot.com/2009/01/viterbi-algorithm-in-c-and-using-stl.html
   by using BOOST and cleaning-up the code
 */


#include <iostream>
#include <string>
#include <vector>
#include <boost/unordered_map.hpp>


namespace Viterbi
{
  template< class TState, class TProbability >
  class HMM
  {
  public:

    typedef HMM< TState, TProbability > Self;

    typedef TState                    StateType;
    typedef std::vector< StateType >  StateVectorType;

    typedef TProbability              ProbabilityType;

    typedef boost::unordered_map< StateType, ProbabilityType>       StateProbabilityMapType;
    typedef boost::unordered_map< StateType, StateProbabilityMapType> StateStateProbabilityMapType;

    HMM(){}
    ~HMM(){}

    void SetStates( const StateVectorType& states )
    { m_States = states; }

    const StateVectorType& GetStates() const
    { return m_States; }


    void SetObservations( const StateVectorType& observations )
    { m_Observations = observations; }

    const StateVectorType& GetObservations() const
    { return m_Observations;}


    void SetStartProbability( const StateProbabilityMapType& prob )
    { m_StartProbability = prob; }

    StateProbabilityMapType& GetStartProbability()
    { return m_StartProbability; }


    void SetTransitionProbability( const StateStateProbabilityMapType& t )
    { m_TransitionProbability = t; }

    StateStateProbabilityMapType& GetTransitionProbability()
    { return m_TransitionProbability; }


    void SetEmissionProbability( const StateStateProbabilityMapType& e )
    { m_EmissionProbability = e; }

    StateStateProbabilityMapType& GetEmissionProbability()
    { return m_EmissionProbability; }

  private:
    HMM( const HMM& );
    void operator = ( const HMM& );

    StateVectorType               m_States;
    StateVectorType               m_Observations;
    StateProbabilityMapType       m_StartProbability;
    StateStateProbabilityMapType  m_TransitionProbability;
    StateStateProbabilityMapType  m_EmissionProbability;
  };

  template< class TState, class TProbability >
  std::ostream& operator << (std::ostream& os, HMM< TState, TProbability > &);

  //
  // computes total probability for observation
  // most likely viterbi path
  // and probability of such path
  //
  template< class TState, class TProbability >
  class ForwardViberti
  {
  public:

    typedef TState                    StateType;
    typedef std::vector< StateType >  StateVectorType;

    typedef TProbability              ProbabilityType;

    typedef boost::unordered_map< StateType, ProbabilityType>         StateProbabilityMapType;
    typedef boost::unordered_map< StateType, StateProbabilityMapType> StateStateProbabilityMapType;

    typedef HMM< TState, TProbability > HMMType;

    ForwardViberti( HMMType& hmm )
    {
      this->setHMM( hmm );

      m_TotalProbabilityOfObservationSequence = 0;
      m_TotalProbabilityOfVibertiPath         = 0;
    }

    ~ForwardViberti() {}

    void Update()
    {
      forward_viterbi();
    }

    ProbabilityType GetTotalProbabilityOfObservationSequence() const
    { return m_TotalProbabilityOfObservationSequence; }

    ProbabilityType GetTotalProbabilityOfVibertiPath() const
    { return m_TotalProbabilityOfVibertiPath; }

    const StateVectorType& GetFinalPath() const
    { return m_FinalPath; }

  private:
    ForwardViberti( const ForwardViberti& );
    void operator = ( const ForwardViberti& );

    StateVectorType               m_Observations;
    StateVectorType               m_States;
    StateProbabilityMapType       m_StartProbability;
    StateStateProbabilityMapType  m_TransitionProbability;
    StateStateProbabilityMapType  m_EmissionProbability;

    ProbabilityType               m_TotalProbabilityOfObservationSequence;
    ProbabilityType               m_TotalProbabilityOfVibertiPath;
    StateVectorType               m_FinalPath;


    void setHMM( HMMType& hmm )
    {
      m_Observations          = hmm.GetObservations();
      m_States                = hmm.GetStates();
      m_StartProbability      = hmm.GetStartProbability();
      m_TransitionProbability = hmm.GetTransitionProbability();
      m_EmissionProbability   = hmm.GetEmissionProbability();
    }

    void forward_viterbi();
  };
}

#include "viterbi.hpp"
#endif
