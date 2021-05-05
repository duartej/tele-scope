
// Daniel Pitzl, Mar 2020
// thr scan: eff

// root -l thr.C

//void thr()
{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.0, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.22 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

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

  //--------------------------------------------------------------------
  // square canvas:
  //            topleft x, y, width x, y
  TCanvas c1("c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "thr.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  string HASH {"#"};

  //----------------------------------------------------------------------------
  // chip 577 FDB150P RD53A 100x25 P1 default design, KaZ 1E16 n_eq/cm2 on Zh card
  // thr 1160

  cout << "try to open dat file for c577";
  ifstream Dstream77( "thr577.dat" );
  if( !Dstream77 ) {
    cout << ": failed" << endl;
    return;
  }
  cout << ": succeed" << endl;

  // Read file by lines:

  vector <double> vt77;
  vector <double> ef77;
  vector <double> vh77;

  while( Dstream77.good() && ! Dstream77.eof() ) {

    string rl;
    getline( Dstream77, rl ); // read one line  = event into string
    if( rl.empty() ) continue;
    if( rl.substr(0,1) == HASH ) // comments start with #
      continue;

    istringstream strm( rl ); // tokenize string

    double thr;
    int run;
    double eff;
    int hot;

    strm >> thr;
    strm >> run;
    strm >> eff;
    strm >> hot;

    vt77.push_back(thr);
    ef77.push_back(eff);
    vh77.push_back(hot);

  } // while lines

  int n77 = vt77.size();
  cout << "    " << n77 << " runs for chip 577i" << endl;

  //----------------------------------------------------------------------------
  // eff vs thr

  heff = new
    TH1F( "heff",
	  "RD53A;threshold [e];RD53A Lin efficiency [%]",
	  29, 0, 2900 ); // axis range
  heff->SetStats(kFALSE); // no statistics
  heff->SetMinimum(80);
  heff->SetMaximum(100);
  heff->Draw();

  geff = new TGraph( n77, &vt77[0], &ef77[0] );
  geff->SetLineColor(4);
  geff->SetMarkerColor(4);
  geff->SetMarkerStyle(20);
  geff->SetMarkerSize(1.4);
  geff->Draw("P"); // without axis option: overlay

  TLegend * leff = new TLegend( 0.2, 0.2, 0.7, 0.28 );
  leff->AddEntry( geff, "577: 9.9#upoint10^{15} n_{Ka}, 400 V ", "p" );
  leff->Draw( "same" );

  c1.Print( "thr.pdf" );

  //----------------------------------------------------------------------------
  // hot vs thr

  hhot = new
    TH1F( "hhot",
	  "RD53A;threshold [e];RD53A Lin hot pixels",
	  29, 0, 2900 ); // axis range
  hhot->SetStats(kFALSE); // no statistics
  gPad->SetLogy(1);
  hhot->SetMinimum(1);
  hhot->SetMaximum(1010);
  hhot->Draw();

  ghot = new TGraph( n77, &vt77[0], &vh77[0] );
  ghot->SetLineColor(2);
  ghot->SetMarkerColor(2);
  ghot->SetMarkerStyle(29);
  ghot->SetMarkerSize(1.9);
  ghot->Draw("P"); // without axis option: overlay

  TLegend * lhot = new TLegend( 0.2, 0.2, 0.7, 0.28 );
  lhot->AddEntry( ghot, "577: 9.9#upoint10^{15} n_{Ka}, 400 V ", "p" );
  lhot->Draw( "same" );

  c1.Print( "thr.pdf" );
  gPad->SetLogy(0);

  //----------------------------------------------------------------------------
  // done:

  c1.Print( "thr.pdf]" ); // ] closes file
  cout << "evince thr.pdf" << endl;

}
