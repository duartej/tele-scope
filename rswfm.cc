
// Daniel Pitzl, DESY, Sep 2019
// read Rohde & Schwarz waveform file

#include <iostream> // cout
#include <fstream> // files

#include <TFile.h> // ROOT
#include <TH1.h>
#include <TH2.h>
#include <TProfile.h>

using namespace std;

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "usage: rswfm run_0.Wfm.bin" << endl;
    return 1;
  }

  // file = last arg

  string fileName( argv[argc-1] );

  cout << "try to open: " << fileName;

  ifstream wff( fileName, ios::binary );

  if (! wff ) {
    cout << " failed" << endl;
    return 1;
  }
  cout << " : success" << endl << flush;

  TFile * histoFile = new TFile( "rswfm.root", "RECREATE" );

  // book histos:

  TH1I hdtus( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms( "dtms", "time between events;time between events [ms];events", 100, 0, 100 );
  TH1I t1( "t1", "events vs time;time [s];events/bin", 100, 0, 1 );
  TH1I t2( "t2", "events vs time;time [s];events/bin", 10*12.5, 0, 10 );

  TProfile wave( "wave", "average wave forms;sample;<ADC>", 8000, -0.5, 7999.5 );
  TH1I hadc( "adc", "ADC;ADC;samples", 256, -128.5, 127.5 );
  TH1I hpeak( "peak", "peak ADC;peak ADC;waves", 256, -128.5, 127.5 );
  TH1I hpeakpos( "peakpos", "any peak position;peak sample index;waves", 8000, -0.5, 7999.5 );
  TH1I hpulsepos( "pulsepos", "pulse peak position;peak sample index;peak waves", 8000, -0.5, 7999.5 );
  TProfile pulse( "pulse", "average pulse;sample;<ADC>", 8000, -0.5, 7999.5 );

  //----------------------------------------------------------------------------

  bool ldb = 0; // debug flag all
  bool ldb2 = 0; // debug flag times

  int nby = 0; // bytes
  int nev = 0; // triggers

  double ts0 = 0;
  double prevts = 0;

  //----------------------------------------------------------------------------

  char y;
  wff.read( ( char * ) &y, sizeof( y ) ); // read aheaad

  //while( !wff.eof() && nby <= 80999 ) {
  while( !wff.eof() ) { 

    /*
      one trigger:

         8 bytes pre
         8 bytes timestamp
         7 bytes pre
      8000 bytes waveform
         17 bytes post
      ----
      8040
    */

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // pre-bytes:

    if( ldb ) cout << wff.tellg();

    for( int k = 0; k < 7; ++k ) {
      wff.read( ( char * ) &y, sizeof( y ) );
      if( ldb ) cout << " " << (int) y;
      ++nby;
    }
    if( ldb ) cout << "  " << nby << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // time stamp:

    double ts;
    wff.read( ( char * ) &ts, sizeof( ts ) );

    if( ldb || ldb2 ) cout << ts << " s" << endl;

    nby += sizeof( ts );

    if( nev == 0 )
      ts0 = ts; // most negative

    t1.Fill( ts-ts0 );
    t2.Fill( ts-ts0 );

    double dt = ts-prevts;
    if( nev > 0 ) {
      hdtus.Fill( dt * 1E6 ); // [us]
      hdtms.Fill( dt * 1E3 ); // [ms]
    }
    prevts = ts; // remember

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // more pre-bytes:

    if( ldb ) cout << wff.tellg();

    for( int k = 0; k < 7; ++k ) {

      wff.read( ( char * ) &y, sizeof( y ) );
      if( ldb ) cout << " " << (int) y;
      ++nby;
    }
    if( ldb ) cout << "  " << nby << endl;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // waveform:

    if( ldb ) cout << wff.tellg();

    int ymin = 128;
    int kmin = 0;
    char yy[8000];

    for( int k = 0; k < 8000; ++k ) { // from xml RecordLength

      wff.read( ( char * ) &y, sizeof( y ) );
      ++nby;
      yy[k] = y;
      hadc.Fill( y ); // mean 25.66, rms 2.4
      wave.Fill( k, y );
      if( y < ymin ) {
	ymin = y;
	kmin = k;
      }
    } // k samples
    if( ldb ) cout << "  " << nby << endl;

    hpeak.Fill( ymin );
    hpeakpos.Fill( kmin );

    if( ymin < 0 ) { // pulse

      hpulsepos.Fill( kmin );

      for( int k = 0; k < 8000; ++k ) { // from xml RecordLength
	pulse.Fill( k, yy[k] );
      }

    } // peak

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // post bytes:

    if( ldb ) cout << wff.tellg();
    for( int k = 0; k < 17; ++k ) {
      wff.read( ( char * ) &y, sizeof( y ) );
      if( ldb ) cout << " " << (int) y;
      ++nby;
    }
    if( ldb ) cout << "  " << nby << endl;

    ++nev;
    if( nev%1000 == 0 ) cout << " " << nev << flush;

    wff.read( ( char * ) &y, sizeof( y ) ); // read ahead, or EOF

  } // file

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  cout << endl;
  cout << "events " << nev << endl
       << "bytes  " << nby << endl
       << "time0  " << ts0 << endl;

  wff.close(); 

  cout << endl << histoFile->GetName() << endl;
  cout << endl;
  histoFile->ls();

  histoFile->Write();
  histoFile->Close();

  cout << endl;

  return 0;
}
