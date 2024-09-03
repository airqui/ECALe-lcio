#include "MergeEvents.h"
#include <iostream>

#include <marlin/Global.h>
#include "marlin/ProcessorEventSeeder.h"

#ifdef MARLIN_USE_AIDA
#include <marlin/AIDAProcessor.h>

#include <AIDA/IHistogramFactory.h>
#include <AIDA/ICloud1D.h>
//#include <AIDA/IHistogram1D.h>
#endif

#include "EVENT/LCEvent.h"
#include "IO/LCReader.h"
#include "IO/LCWriter.h"
#include <EVENT/LCCollection.h>
#include <IMPL/LCCollectionVec.h>
#include <EVENT/MCParticle.h>
#include "IO/LCReader.h"
#include "UTIL/LCTOOLS.h"
#include "Merger.h"

#include "CLHEP/Random/RandPoisson.h"
#include "CLHEP/Random/RandFlat.h"
// #include <time.h>

using namespace lcio ;
using namespace marlin ;

namespace MergeEvents {
  
  MergeEvents aMergeEvents ;
  
  /// Set the lcio file name
  void LCFileHandler::setFileName(const std::string& fname) {
    _fileName = fname ;
  }
  
  //===========================================================================================================================
  
  /// Get the number of events available in the file
  unsigned int LCFileHandler::getNumberOfEvents() const {
    openFile() ;
    return _lcReader->getNumberOfEvents() ;
  }
  
  //===========================================================================================================================
  
  /// Get the event number at the specified index (look in the event map)
  unsigned int LCFileHandler::getEventNumber(unsigned int index) const {
    openFile() ;
    return _eventMap.at( index * 2 + 1 ) ;
  }
  
  //===========================================================================================================================
  
  /// Get the run number at the specified index (look in the event map)
  unsigned int LCFileHandler::getRunNumber(unsigned int index) const {
    openFile() ;
    return _eventMap.at( index * 2 ) ;
  }
  
  //===========================================================================================================================
  
  /// Read the specified event, by run and event number
  EVENT::LCEvent* LCFileHandler::readEvent(int runNumber, int eventNumber) {
    openFile();
    streamlog_out( DEBUG6 ) << "*** Reading event from file : '" << _fileName 
          << "',  event number " << eventNumber << " of run " << runNumber << "." << std::endl ;
    return _lcReader->readEvent( runNumber, eventNumber, LCIO::UPDATE );
  }
  
  //===========================================================================================================================
  
  /// Proxy method to open the LCIO file
  void LCFileHandler::openFile() const {
    if(nullptr == _lcReader) {
      _lcReader = std::shared_ptr<IO::LCReader>( LCFactory::getInstance()->createLCReader( LCReader::directAccess ) ) ;
      
      streamlog_out( MESSAGE ) << "*** Opening file for MergeEvents, file name:" << _fileName << std::endl ;
            
      _lcReader->open( _fileName ) ;
      _lcReader->getEvents( _eventMap ) ;
      
      streamlog_out( MESSAGE ) << "*** Opening file for MergeEvents : number of available events: " << _lcReader->getNumberOfEvents() << std::endl ;
    }
  }

  //===========================================================================================================================
  //===========================================================================================================================
  
  MergeEvents::MergeEvents() : Processor("MergeEvents") {
    
    // modify processor description
    _description = "Opens a second (chain of) lcio file(s) and MergeEventss events..." ;
  
    // register steering parameters: name, description, class-variable, default value
    registerProcessorParameter( "InputFileNames" , 
				"Name of the lcio input file(s)"  ,
				_fileNames ,
				StringVec( {"undefined.slcio"} ) ) ;
  
    /*registerProcessorParameter( "NumberMergeEvents" , 
				"MergeEvents each event with this number of background events. (default 0)" ,
				_numMergeEvents ,
				static_cast<int>(0) ) ; */
     
    registerProcessorParameter( "CollectionMap" , 
				"Pairs of collection to be merged"  ,
				_MergeEventsCollections ,
				StringVec( {"MCParticle", "MCParticle"} ) ) ;
        
    registerProcessorParameter( "ExcludeCollections" , 
        "List of collections to exclude for merging"  ,
        _excludeCollections ,
        StringVec() ) ;
  }
  
  //===========================================================================================================================
  
  const std::string & MergeEvents::name() const {
    return Processor::name() ;
  }

  //===========================================================================================================================

  void MergeEvents::init() { 
    
    if( ! LCIO_VERSION_GE( 1 , 4 ) ) {
      throw Exception("  MergeEvents requires LCIO v1.4 or higher \n"
		      "  - please upgrade your LCIO version or disable the MergeEvents processor ! ") ;
    }

    // usually a good idea to
    printParameters() ;
    
    // prepare the lcio file handlers
    _lcFileHandlerList.resize( _fileNames.size() ) ;
    
    for ( unsigned int i=0 ; i<_fileNames.size() ; i++ )
      _lcFileHandlerList.at( i ).setFileName( _fileNames.at( i ) ) ;
  
    // initalisation of random number generator
    Global::EVENTSEEDER->registerProcessor(this) ;

    // preparing collection map for merge
    StringVec::iterator endIter = _MergeEventsCollections.end() ;
  
    if ( _MergeEventsCollections.size() % 2 ) { 
      streamlog_out( WARNING ) << "Odd number of collection names, last collection ignored." << std::endl ;
      --endIter ;
    }
  
    // treating pairs of collection names
    for ( auto iter =_MergeEventsCollections.begin() ; iter != endIter ; ++iter ) {  
      std::string key = (*iter) ;
      ++iter ;
      _MergeEventsCollectionMap[key] = (*iter) ;
    }

    _nRun = 0 ;
    _nEvt = 0 ;  
  }

  //===========================================================================================================================

  void MergeEvents::processRunHeader( LCRunHeader* ) {
  
    _nRun++ ;
  } 

  //===========================================================================================================================

  void MergeEvents::modifyEvent( LCEvent * evt ) {

    if( isFirstEvent() ) {
      // get it here and not in init as files are opened on function call
      _nAvailableEvents = getNAvailableEvents() ;
      
      streamlog_out( MESSAGE ) << "MergeEvents::modifyEvent: total number of available events to MergeEvents: " << _nAvailableEvents << std::endl ;
    }

    // initalisation of random number generator
    int eventSeed = Global::EVENTSEEDER->getSeed(this);
    CLHEP::HepRandom::setTheSeed( eventSeed );


    //unsigned int nEventsToMergeEvents = _numMergeEvents ;
  
   /* if ( parameterSet("expBG") ) {
      nEventsToMergeEvents += CLHEP::RandPoisson::shoot( _expBG ) ;
    }*/

    streamlog_out( DEBUG6 ) << "** Processing event nr " << evt->getEventNumber() << " run " <<  evt->getRunNumber() 
     			//    << "\n**  Merge initial plus " << nEventsToMergeEvents << " added events." 
			    //<< " ( seeded CLHEP::HepRandom with seed = " << eventSeed  << ") " 
			    << std::endl;
  
    int nOverlaidEvents(0);
    EVENT::FloatVec overlaidEventIDs, overlaidRunIDs;
    
    //for(unsigned int i=0 ; i < nEventsToMergeEvents ; i++ ) {

      EVENT::LCEvent *MergeEventsEvent = readNextEvent(_nEvt) ;

      if( nullptr != MergeEventsEvent ) {
	       
      
      overlaidEventIDs.push_back( MergeEventsEvent->getEventNumber() );
      overlaidRunIDs.push_back( MergeEventsEvent->getRunNumber() );
      
      ++nOverlaidEvents ;

      streamlog_out( DEBUG6 ) << " will MergeEvents event " << MergeEventsEvent->getEventNumber() << " - run " << MergeEventsEvent->getRunNumber() << std::endl ;

      std::map<std::string, std::string> collectionMap;
      
      if ( ( _MergeEventsCollectionMap.empty() || ! parameterSet("CollectionMap") ) ) {
        auto collectionNames = MergeEventsEvent->getCollectionNames();
        for ( auto collection : *collectionNames ) {
          collectionMap[collection] = collection;
          streamlog_out( DEBUG6 ) << "Collection map -> " << collection << std::endl;
        }
      }
      else {
        collectionMap = _MergeEventsCollectionMap;
      }
      
      // Remove collections to exclude from the collection map
      if(not _excludeCollections.empty()) {
        for ( auto excludeCol : _excludeCollections ) {
          auto findIter = collectionMap.find( excludeCol );
          if( collectionMap.end() != findIter ) {
            collectionMap.erase( findIter );
          }
        }
      }

	     Merger::merge( MergeEventsEvent, evt, &collectionMap );
    //}
    
    _nTotalMergeEventsEvents += nOverlaidEvents;
    
    // Write info to event parameters
    std::string paramName = "MergeEvents." + this->name() + ".nEvents";
    evt->parameters().setValue(paramName, nOverlaidEvents);
    paramName = "MergeEvents." + this->name() + ".eventIDs";
    evt->parameters().setValues(paramName, overlaidEventIDs);
    paramName = "MergeEvents." + this->name() + ".runIDs";
    evt->parameters().setValues(paramName, overlaidRunIDs);
    
    int totalMergeEvents = evt->parameters().getIntVal("MergeEvents.nTotalEvents"); // returns 0 if the key doesn't exists
    totalMergeEvents += nOverlaidEvents;
    evt->parameters().setValue("MergeEvents.nTotalEvents", totalMergeEvents);

    _nEvt ++ ;
      } else {
        streamlog_out( ERROR ) << _nEvt<<" ++++++++++ Nothing to MergeEvents +++++++++++ \n " ;
	      //_nEvt ++ ;
        //_nTotalMergeEventsEvents += nOverlaidEvents;
      } 
  }

  //===========================================================================================================================

  void MergeEvents::end() { 

    streamlog_out( MESSAGE ) << " ------------------------------------------ " 
			     << "   MergeEvents processor " << _nTotalMergeEventsEvents << " events on " << _nEvt << " physics events."
			     << std::endl ;
  }


  //===========================================================================================================================

  EVENT::LCEvent* MergeEvents::readNextEvent(int eventIndex) {
    
    EVENT::LCEvent* MergeEventsEvent = nullptr ;
    
    // get the event index to random pick an event among the possible files, in case the second file is smaller than the first, we merge random events

    unsigned int currentEventIndex(0);
    
    streamlog_out( DEBUG ) << "MergeEvents::readNextEvent: index = " << eventIndex  << " over " << _nAvailableEvents << std::endl ;
    
    for ( auto &handler : _lcFileHandlerList ) {

      std::cout<<"e o"<<std::endl;
      if(eventIndex >= _nAvailableEvents) {
        const unsigned int eventIndex_random = CLHEP::RandFlat::shoot( static_cast<double>( _nAvailableEvents  ) ) ;
        streamlog_out( DEBUG )<< "MergeEvents::readNextEvent: the second file has size:"<< _nAvailableEvents << 
        " which is smaller than the current eventN of the first file: " << eventIndex << 
        "  so we merge a random event --> evt: " <<  eventIndex_random << std::endl ;
        eventIndex=eventIndex_random;
      }

      if ( currentEventIndex <= eventIndex && eventIndex <= currentEventIndex +  _nAvailableEvents ) {
        
        const int eventNumber = handler.getEventNumber( eventIndex);//-currentEventIndex ) ;
        const int runNumber = handler.getRunNumber( eventIndex);//-currentEventIndex ) ;
        
        MergeEventsEvent = handler.readEvent( runNumber, eventNumber ) ;
        
        if( nullptr == MergeEventsEvent ) {
          streamlog_out( ERROR ) << "MergeEvents::readNextEvent: Could not read event " << eventNumber << "  from  run " <<  runNumber << std::endl ;
        }
        break;
      }
      currentEventIndex += handler.getNumberOfEvents() ;
    }
    
    return MergeEventsEvent ;
  }
  
  //===========================================================================================================================
  
  unsigned int MergeEvents::getNAvailableEvents() const
  {
    unsigned int totalNEvents(0);
    
    for ( auto &handler : _lcFileHandlerList ) {
      totalNEvents += handler.getNumberOfEvents();
    }
    
    return totalNEvents;
  }

} // end namespace 
