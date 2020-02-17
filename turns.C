
// Daniel Pitzl, May 2019
// RFD53A turn scan fresh and irrad
// root -l turns.C
/*
#include <TCanvas.h>
#include <TGraph.h>
#include <TH1.h>
#include <TF1.h>

//------------------------------------------------------------------------------
void turns()*/
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleOffset( 1.5, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(111111);
  gStyle->SetStatFormat("8.6g"); // more digits, default is 6.4g
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize(1); // no 'shadow'

  gStyle->SetStatX(0.95);
  gStyle->SetStatY(0.90);

  gStyle->SetPalette(1); // rainbow colors

  gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //                topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );

  c1.Print( "turns.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  string HASH( "#" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read fresh data:

  cout << "try to open turn524.dat";
  ifstream D524( "turn524.dat" );
  if( !D524 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vt24;
  vector <double> vq24;
  vector <double> vn24;
  vector <double> vr24;

  while( D524.good() && ! D524.eof() ) {

    string rl;
    getline( D524, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int run;
    double trn;
    double qln;
    double qdf;
    double npx;
    double dxc;

    strm >> run;
    strm >> trn;
    strm >> qln;
    strm >> qdf;
    strm >> npx;
    strm >> dxc;

    vt24.push_back(trn);
    vq24.push_back(qln);
    vn24.push_back(npx);
    vr24.push_back(dxc);

  } // while lines

  int n524 = vt24.size();
  cout << n524 << " turns for chip 524" << endl;

  vector <double> vs24(n524);
  for( int i = 0; i < n524; ++i ) {
    vs24[i] = sqrt( vr24[i]*vr24[i] - 6.5*6.5 ); // subtract track sixdxc/2
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read PS data:

  cout << "try to open turn511.dat";
  ifstream D511( "turn511.dat" );
  if( !D511 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vt11;
  vector <double> ve11;
  vector <double> vn11;
  vector <double> vq11;

  while( D511.good() && ! D511.eof() ) {

    string rl;
    getline( D511, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int run;
    double trn;
    double eff;
    double col;
    double q;

    strm >> run;
    strm >> trn;
    strm >> eff;
    strm >> col;
    strm >> q;

    vt11.push_back(trn);
    ve11.push_back(eff*0.01);
    vn11.push_back(col);
    vq11.push_back(q);

  } // while lines

  int n511 = vt11.size();
  cout << n511 << " turns for chip 511" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // read KaZ data:

  cout << "try to open turn521.dat";
  ifstream D521( "turn521.dat" );
  if( !D521 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vt21;
  vector <double> ve21;
  vector <double> vn21;
  vector <double> vq21;
  vector <double> sx21;
  vector <double> sy21;

  while( D521.good() && ! D521.eof() ) {

    string rl;
    getline( D521, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int run;
    double trn;
    double eff;
    double col;
    double q;
    double dxc;
    double dyc;

    strm >> run;
    strm >> trn;
    strm >> eff;
    strm >> col;
    strm >> q;
    strm >> dxc;
    strm >> dyc;

    vt21.push_back(trn);
    ve21.push_back(eff*0.01);
    vn21.push_back(col);
    vq21.push_back(q);
    sx21.push_back( sqrt( dxc*dxc - 7.86*7.86 ) ); // subtract track
    sy21.push_back( sqrt( dyc*dyc - 7.86*7.86 ) ); // subtract track

  } // while lines

  int n521 = vt21.size();
  cout << n521 << " turns for chip 521" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // 547i:

  cout << "try to open turn547.dat";
  ifstream D547( "turn547.dat" );
  if( !D547 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vt47; // turn
  vector <double> rx47; // res
  vector <double> sx47; // sig
  vector <double> vn47; // ncol

  while( D547.good() && ! D547.eof() ) {

    string rl;
    getline( D547, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    int run;
    double trn;
    double dxc;
    double col;

    strm >> run;
    strm >> trn;
    strm >> dxc;
    strm >> col;

    vt47.push_back(-trn); // flip sign
    vn47.push_back(col);
    rx47.push_back( dxc ); // residual width
    double tel = 7.3; // telescope
    if( run >= 38111 ) // larger opening
      //tel = 9.3; // [mu]
      tel = 10.4; // [mu]
    sx47.push_back( sqrt( dxc*dxc - tel*tel ) ); // subtract track

  } // while lines

  int n547 = vt47.size();
  cout << n547 << " turns for chip 547" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // eff vs turn:

  TH1F * he = new
    TH1F( "he",
	  "efficiency vs turn;turn angle [deg];Lin efficiency",
	  31, 0, 31 ); // axis range
  he->SetStats(kFALSE); // no statistics
  he->SetMinimum(0.9);
  he->SetMaximum(1);
  he->Draw();

  TGraph * ge21 = new TGraph( n521, &vt21[0], &ve21[0] );
  ge21->SetMarkerColor(kGreen+3);
  ge21->SetMarkerStyle(21);
  ge21->SetMarkerSize(1.5);
  ge21->Draw("P"); // without axis option: overlay

  TGraph * ge11 = new TGraph( n511, &vt11[0], &ve11[0] );
  ge11->SetMarkerColor(1);
  ge11->SetMarkerStyle(22);
  ge11->SetMarkerSize(1.5);
  ge11->Draw("P"); // without axis option: overlay

  TLegend * lge = new TLegend( 0.36, 0.20, 0.93, 0.32 );
  lge->AddEntry( ge21, "521 5.3#upoint10^{15} n_{Ka}/cm^{2}, 775 V ", "p" );
  lge->AddEntry( ge11, "511 6.6#upoint10^{15} n_{PS}/cm^{2}, 750 V ", "p" );
  lge->Draw( "same" );

  c1.Print( "turns.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ncol vs turn:

  gStyle->SetTitleOffset( 1.7, "y" );

  TH1F * hn = new
    TH1F( "hn",
	  "cluster size vs turn;turn angle [deg];Lin cluster size [columns]",
	  31, 0, 31 ); // axis range
  hn->SetStats(kFALSE); // no statistics
  hn->SetMinimum(1);
  hn->SetMaximum(3);
  hn->Draw();

  TGraph * gn24 = new TGraph( n524, &vt24[0], &vn24[0] );
  gn24->SetMarkerColor(4);
  gn24->SetMarkerStyle(20);
  gn24->SetMarkerSize(1.5);
  gn24->Draw("P"); // without axis option: overlay

  TGraph * gn21 = new TGraph( n521, &vt21[0], &vn21[0] );
  gn21->SetMarkerColor(kGreen+3);
  gn21->SetMarkerStyle(21);
  gn21->SetMarkerSize(1.5);
  gn21->Draw("P"); // without axis option: overlay

  TGraph * gn11 = new TGraph( n511, &vt11[0], &vn11[0] );
  gn11->SetMarkerColor(1);
  gn11->SetMarkerStyle(22);
  gn11->SetMarkerSize(1.5);
  gn11->Draw("P"); // without axis option: overlay

  TLegend * lgn = new TLegend( 0.17, 0.74, 0.72, 0.88 );
  lgn->AddEntry( gn24, "524 fresh, 120 V", "p" );
  lgn->AddEntry( gn21, "521 5.3#upoint10^{15} n_{Ka}/cm^{2}, 775 V ", "p" );
  lgn->AddEntry( gn11, "511 6.6#upoint10^{15} n_{PS}/cm^{2}, 750 V ", "p" );
  lgn->Draw( "same" );

  c1.Print( "turns.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // q vs turn:

  TH1F * hq = new
    TH1F( "hq",
	  "cluster ToT vs turn;turn angle [deg];Lin <cluster charge> [ToT]",
	  31, 0, 31 ); // axis range
  hq->SetStats(kFALSE); // no statistics
  hq->SetMinimum( 0);
  hq->SetMaximum(20);
  hq->Draw();

  // fresh:

  TGraph * gq24 = new TGraph( n524, &vt24[0], &vq24[0] );
  gq24->SetMarkerColor(4);
  gq24->SetMarkerStyle(20);
  gq24->SetMarkerSize(1.5);
  gq24->Draw("P"); // without axis option: overlay

  // irrad:

  TGraph * gq21 = new TGraph( n521, &vt21[0], &vq21[0] );
  gq21->SetMarkerColor(kGreen+3);
  gq21->SetMarkerStyle(21);
  gq21->SetMarkerSize(1.5);
  gq21->Draw("P"); // without axis option: overlay

  TGraph * gq11 = new TGraph( n511, &vt11[0], &vq11[0] );
  gq11->SetMarkerColor(1);
  gq11->SetMarkerStyle(22);
  gq11->SetMarkerSize(1.5);
  gq11->Draw("P"); // without axis option: overlay

  TLegend * lgq = new TLegend( 0.36, 0.17, 0.93, 0.31 );
  lgq->AddEntry( gq24, "524 fresh, 120 V", "p" );
  lgq->AddEntry( gq21, "521 5.3#upoint10^{15} n_{Ka}/cm^{2}, 775 V ", "p" );
  lgq->AddEntry( gq11, "511 6.6#upoint10^{15} n_{PS}/cm^{2}, 750 V ", "p" );
  lgq->Draw( "same" );

  c1.Print( "turns.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // sig vs turn

  TH1F * hs = new
    TH1F( "hs",
	  "resolution vs turn;turn angle [deg];Lin resolution [#mum]",
	  49, -10, 39 ); // axis range
  hs->SetStats(kFALSE); // no statistics
  hs->SetMinimum(0);
  hs->SetMaximum(16);
  hs->Draw();

  TGraph * gs24 = new TGraph( n524, &vt24[0], &vs24[0] );
  gs24->SetLineColor(4);
  gs24->SetMarkerColor(4);
  gs24->SetMarkerStyle(20);
  gs24->SetMarkerSize(1.5);
  gs24->Draw("Pc"); // without axis option: overlay

  TGraph * gs47 = new TGraph( n547, &vt47[0], &sx47[0] );
  gs47->SetLineColor(1);
  gs47->SetMarkerColor(1);
  gs47->SetMarkerStyle(21);
  gs47->SetMarkerSize(1.5);
  gs47->Draw("Pc"); // without axis option: overlay

  TGraph * gs21 = new TGraph( n521, &vt21[0], &sx21[0] );
  gs21->SetLineColor(kGreen+3);
  gs21->SetMarkerColor(kGreen+3);
  gs21->SetMarkerStyle(21);
  gs21->SetMarkerSize(1.5);
  gs21->Draw("Pc"); // without axis option: overlay

  TLegend * lgs = new TLegend( 0.36, 0.16, 0.93, 0.32 );
  lgs->AddEntry( gs24, "524 fresh, 120 V", "p" );
  lgs->AddEntry( gs47, "547 4.4#upoint10^{15} n_{Ka}/cm^{2}, 764 V ", "p" );
  lgs->AddEntry( gs21, "521 5.3#upoint10^{15} n_{Ka}/cm^{2}, 775 V ", "p" );
  lgs->Draw( "same" );

  c1.Print( "turns.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "turns.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf turns.ps" );
  ierr = system( "rm -f  turns.ps" );
  cout << "evince turns.pdf" << endl;

}
