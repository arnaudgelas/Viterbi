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

class MyException : public std::exception
{
public:
  MyException( const std::string& msg ) : Message( msg ) {}
  MyException( const char* msg ) : Message( msg ) {}
  MyException( const MyException& e ) : Message( e.Message ) {}
  virtual ~MyException() throw() {}

  virtual const char *what() const throw()
  {
    return Message.c_str();
  }

private:
  void operator = ( const MyException& );

  std::string Message;
};

namespace Viterbi
{
  template< class TState, class TObservation, class TProbability >
  class HMM
  {
  public:

    typedef HMM< TState, TObservation, TProbability > Self;

    typedef TState                    StateType;
    typedef std::vector< StateType >  StateVectorType;

    typedef TObservation                    ObservationType;
    typedef std::vector< ObservationType >  ObservationVectorType;

    typedef TProbability              ProbabilityType;

    typedef boost::unordered_map< StateType, ProbabilityType>               StateProbabilityMapType;
    typedef boost::unordered_map< ObservationType, ProbabilityType>         ObservationProbabilityMapType;

    typedef boost::unordered_map< StateType, StateProbabilityMapType>       StateStateProbabilityMapType;
    typedef boost::unordered_map< StateType, ObservationProbabilityMapType> StateObservationProbabilityMapType;

    HMM(){}
    ~HMM(){}

    void SetStates( const StateVectorType& states )
    { m_States = states; }

    const StateVectorType& GetStates() const
    { return m_States; }


    void SetObservations( const ObservationVectorType& observations )
    { m_Observations = observations; }

    const ObservationVectorType& GetObservations() const
    { return m_Observations;}


    void SetStartProbability( const StateProbabilityMapType& prob )
    { m_StartProbability = prob; }

    StateProbabilityMapType& GetStartProbability()
    { return m_StartProbability; }


    void SetTransitionProbability( const StateStateProbabilityMapType& t )
    { m_TransitionProbability = t; }

    StateStateProbabilityMapType& GetTransitionProbability()
    { return m_TransitionProbability; }


    void SetEmissionProbability( const StateObservationProbabilityMapType& e )
    { m_EmissionProbability = e; }

    StateObservationProbabilityMapType& GetEmissionProbability()
    { return m_EmissionProbability; }

  private:
    HMM( const HMM& );
    void operator = ( const HMM& );

    ObservationVectorType               m_Observations;
    StateVectorType                     m_States;
    StateProbabilityMapType             m_StartProbability;
    StateStateProbabilityMapType        m_TransitionProbability;
    StateObservationProbabilityMapType  m_EmissionProbability;
  };

  template< class TState, class TObservation, class TProbability >
  std::ostream& operator << (std::ostream& os, HMM< TState, TObservation, TProbability > &);

  //
  // computes total probability for observation
  // most likely viterbi path
  // and probability of such path
  //
  template< class TState, class TObservation, class TProbability >
  class ForwardViberti
  {
  public:

    typedef TState                    StateType;
    typedef std::vector< StateType >  StateVectorType;

    typedef TObservation                    ObservationType;
    typedef std::vector< ObservationType >  ObservationVectorType;

    typedef TProbability              ProbabilityType;

    typedef boost::unordered_map< StateType, ProbabilityType>               StateProbabilityMapType;
    typedef boost::unordered_map< ObservationType, ProbabilityType>         ObservationProbabilityMapType;

    typedef boost::unordered_map< StateType, StateProbabilityMapType>       StateStateProbabilityMapType;
    typedef boost::unordered_map< StateType, ObservationProbabilityMapType> StateObservationProbabilityMapType;

    typedef HMM< TState, TObservation, TProbability > HMMType;

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

    ObservationVectorType               m_Observations;
    StateVectorType                     m_States;
    StateProbabilityMapType             m_StartProbability;
    StateStateProbabilityMapType        m_TransitionProbability;
    StateObservationProbabilityMapType  m_EmissionProbability;

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
