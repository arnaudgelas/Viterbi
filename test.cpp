#include "viterbi.h"

using namespace Viterbi;

template< class TState, class TObservation, class TProbability >
void init( HMM< TState, TObservation, TProbability >& hmm )
{
  typedef HMM< TState, TObservation, TProbability > HMMType;

  typename HMMType::StateVectorType states;
  states.push_back("Rainy");
  states.push_back("Sunny");
  hmm.SetStates( states );

  typename HMMType::ObservationVectorType observations;
  observations.push_back( "walk");
  observations.push_back( "shop");
  observations.push_back( "clean");
  hmm.SetObservations( observations );

  typename HMMType::StateProbabilityMapType start_probability;
  start_probability["Rainy"] = 0.6;
  start_probability["Sunny"] = 0.4;
  hmm.SetStartProbability( start_probability );

  typename HMMType::StateStateProbabilityMapType transition_probability;
  transition_probability["Rainy"]["Rainy"] = 0.7;
  transition_probability["Rainy"]["Sunny"] = 0.3;
  transition_probability["Sunny"]["Rainy"] = 0.4;
  transition_probability["Sunny"]["Sunny"] = 0.6;
  hmm.SetTransitionProbability( transition_probability );

  typename HMMType::StateObservationProbabilityMapType emission_probability;
  emission_probability["Rainy"]["walk"] = 0.1;
  emission_probability["Rainy"]["shop"] = 0.4;
  emission_probability["Rainy"]["clean"] = 0.5;
  emission_probability["Sunny"]["walk"] = 0.6;
  emission_probability["Sunny"]["shop"] = 0.3;
  emission_probability["Sunny"]["clean"] = 0.1;
  hmm.SetEmissionProbability( emission_probability );
}

int main(int , char* [])
{
  typedef std::string                                         StateType;
  typedef std::string                                         ObservationType;
  typedef double                                              ProbabilityType;
  typedef HMM< StateType, ObservationType, ProbabilityType >  HMMType;

  HMMType h;
  init( h );

  std::cout << h << std::endl;

  ForwardViberti< StateType, ObservationType, ProbabilityType > fv( h );

  try
    {
    fv.Update();
    }
  catch( std::exception& e )
    {
    std::cerr << e.what() << std::endl;
    return EXIT_FAILURE;
    }

  std::cout << std::endl;

  std::cout << "Total probability of the observation sequence: "
      << fv.GetTotalProbabilityOfObservationSequence() << std::endl;

  std::cout << std::endl;

  std::cout << "Probability of the Viterbi path: "
      << fv.GetTotalProbabilityOfVibertiPath() << std::endl;

  std::cout << std::endl;

  std::cout << "The Viterbi path: " << std::endl;

  BOOST_FOREACH( HMMType::StateVectorType::value_type state, fv.GetFinalPath() )
    {
    std::cout << "\tVState: " << state << std::endl;
    }

  return EXIT_SUCCESS;
}
