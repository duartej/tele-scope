
// Daniel Pitzl, DESY, Sep 2019
// read Tektronix waveform file

#include <iostream>
#include <fstream>

using namespace std;

int main()
{
  cout << "hello" << endl;

  ifstream wff( "tek36987.wfm", ios::binary );

  if (! wff ) {
    cout << "Cannot open file";
    return 1;
  }

  char c;
  cout << "length(char) " << sizeof(c) << endl; // 1

  short s;
  cout << "length(shrt) " << sizeof(s) << endl; // 2

  int i;
  cout << "length( int) " << sizeof(i) << endl; // 4

  unsigned u;
  cout << "length(unsigned) " << sizeof(u) << endl; // 4

  float f;
  cout << "length(floa) " << sizeof(f) << endl; // 4

  long l;
  cout << "length(long) " << sizeof(l) << endl; // 8

  double d;
  cout << "length(dble) " << sizeof(d) << endl; // 8

  char vers [8];
  cout << "length(vers) " << sizeof(vers) << endl; // 8

  // file info:

  cout << endl << wff.tellg();
  wff.read( ( char * ) &s, sizeof( s ) );
  cout << ": endian " << s << " = " << hex << s << dec << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &vers, sizeof( vers ) );
  cout << ": version " << vers << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &c, sizeof( c ) );
  cout << ": digits/byte " << c << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": file size " << i << " bytes" << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &c, sizeof( c ) );
  cout << ": bytes/point " << c << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": bytes to curve buffer " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": h zoom " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &f, sizeof( f ) );
  cout << ": h zoom " << f << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &d, sizeof( d ) );
  cout << ": v zoom " << d << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &f, sizeof( f ) );
  cout << ": v zoom " << f << endl;

  cout << wff.tellg();
  char wflab[32];
  wff.read( ( char * ) &wflab, sizeof( wflab ) );
  cout << ": wf label " << wflab << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": last frame " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &s, sizeof( s ) );
  cout << ": size wf header " << s << endl;

  // wave form header:

  cout << endl << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": wf set type " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": wf count " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &l, sizeof( l ) );
  cout << ": counter " << l << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &l, sizeof( l ) );
  cout << ": counter " << l << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": slot ID " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": static flag " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": wf spec cnt " << i << endl;
  int nwf = i;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": dim ref cnt " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": dim ref cnt " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": data type " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &l, sizeof( l ) );
  cout << ": long counter " << l << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": count " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": count " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": curves " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": number of fast frames " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": number of fast frames " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &s, sizeof( s ) );
  cout << ": summary frame " << s << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": px map format " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &l, sizeof( l ) );
  cout << ": px map " << l << endl;

  // Explicit dimension 1 and 2:

  for( int k = 1; k <= 2; ++ k ) {

    cout << endl << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim scale " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim offset " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": dim size " << i << endl;

    cout << wff.tellg();
    char units[20];
    wff.read( ( char * ) &units, sizeof( units ) );
    cout << ": units " << units << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim min " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim max " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim res " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim ref " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": format " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": storage " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": null " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": over " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": under " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": high " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": low  " << i << endl;

    // user:

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": scale " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &units, sizeof( units ) );
    cout << ": units " << units << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": offset " << d << endl;

    cout << wff.tellg();
    //wff.read( ( char * ) &u, sizeof( u ) );
    //cout << ": density " << u << endl;
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": density " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": HRef " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": delay " << d << endl;

  } // dim k

  // Implicit dimensions 1 and 2

  for( int k = 1; k <= 2; ++k ) {

    cout << endl << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim scale " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim offset " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": dim size " << i << endl;

    char units[20];
    cout << wff.tellg();
    wff.read( ( char * ) &units, sizeof( units ) );
    cout << ": units " << units << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim min " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim max " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim res " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": dim ref " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": spacing " << i << endl;

    // user:

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": scale " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &units, sizeof( units ) );
    cout << ": units " << units << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": offset " << d << endl;

    cout << wff.tellg();
    //wff.read( ( char * ) &u, sizeof( u ) );
    //cout << ": density " << u << endl;
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": density " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": HRef " << d << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &d, sizeof( d ) );
    cout << ": delay " << d << endl;

  } // dim k

  // Time Base 1 and 2 information

  for( int k = 1; k <= 2; ++k ) {

    cout << endl << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": spacing " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": sweep " << i << endl;

    cout << wff.tellg();
    wff.read( ( char * ) &i, sizeof( i ) );
    cout << ": base " << i << endl;

  } // dim k

  // Wfm Update specification

  cout << endl;

  cout << endl << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": offset " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &d, sizeof( d ) );
  cout << ": offset " << d << endl;

  double frc;
  cout << wff.tellg();
  wff.read( ( char * ) &frc, sizeof( frc ) );
  int sec;
  wff.read( ( char * ) &sec, sizeof( sec ) );
  cout << ": time " << sec << ". " << frc << endl;

  int sec0 = sec;
  double frc0 = frc;
  int prevsec = sec;
  double prevfrc = frc;

  // curve:

  cout << endl << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": flags " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &i, sizeof( i ) );
  cout << ": check sum type " << i << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &s, sizeof( s ) );
  cout << ": check sum " << s << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &u, sizeof( u ) );
  cout << ": precharge " << u << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &u, sizeof( u ) );
  cout << ": start offset " << u << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &u, sizeof( u ) );
  cout << ": postcharge start " << u << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &u, sizeof( u ) );
  cout << ": postcharge stop  " << u << endl;

  cout << wff.tellg();
  wff.read( ( char * ) &u, sizeof( u ) );
  cout << ": end offset " << u << endl;

  // Wfm Update specification = time stamp

  cout << endl << wff.tellg() << ": time stamps:" << endl;

  for( int n = 1; n < nwf; ++n ) {

    wff.read( ( char * ) &u, sizeof( u ) ); // 4 byte
    if( n==0 ) cout << "offset " << u << endl;

    wff.read( ( char * ) &d, sizeof( d ) );
    if( n==0 ) cout << "offset " << d << endl; // 8 byte

    wff.read( ( char * ) &frc, sizeof( frc ) ); // 8 byte
    wff.read( ( char * ) &sec, sizeof( sec ) ); // 4 byte

    if( n < 99 || n > nwf-99 )
      cout << n << " at " << sec
	   << " + " << frc << " s"
	   << ", dt " << (sec-prevsec + frc-prevfrc)*1e3 << " ms"
	   << endl;
    prevsec = sec;
    prevfrc = frc;
    
  } // n

  // Wfm Curve information

  cout << endl << wff.tellg() << ": Wfm Curve information:" << endl;

  for( int n = 1; n < nwf; ++n ) {

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "flags " << i << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "check sum type " << i << endl;

    wff.read( ( char * ) &s, sizeof( s ) );
    if( n==1 || n==nwf-1 ) cout << "check sum " << s << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "precharge " << i << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "start offset " << i << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "postcharge start " << i << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "postcharge stop  " << i << endl;

    wff.read( ( char * ) &i, sizeof( i ) );
    if( n==1 || n==nwf-1 ) cout << "end offset " << i << endl;

  } // n

  // curve buffer:
  cout << endl << wff.tellg() << ": curve buffer:" << endl;

  int n = 0;
  cout << n;
  while( !wff.eof() && n <= 999 ) { 

    char y;
    wff.read( ( char * ) &y, sizeof( y ) );

    cout << " " << (int) y;

    ++n;
    if( n%10 == 0 ) cout << endl << n;
  }
  cout << endl;
  wff.close(); 

  cout << "run duration was " << sec-sec0 + frc-frc0 << " s"
       << " for " << nwf << " triggers"
       << endl;

  return 0;
}
