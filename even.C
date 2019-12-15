
// Daniel Pitzl, Apr 2019
// RD53A fraction even rowmin
// root -l even.C

//------------------------------------------------------------------------------
//void resvsthr()
{
  // set styles:

  gStyle->SetTextFont(62); // 62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleOffset( 1.4, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

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
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );

  c1.Print( "even.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  gPad->Update();// required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // HLL RD53A 182 Lin:

  const Int_t n = 8;

  // run             35695  35709  35710  35711  35712  35713  35716  35715
  Float_t dac[n] = {   346,   350,   352,   355,   358,   362,   364,   365 }; // threshold DAC
  Float_t thr[n] = {   660,   750,   850,  1000,  1137,  1300,  1376,  1500 }; // threshold e
  Float_t row[n] = { 1.657, 1.564, 1.498, 1.406, 1.348, 1.292, 1.279, 1.255 }; // rows
  Float_t eve[n] = { 68.10, 64.70, 61.60, 57.50, 55.66, 54.00, 53.70, 53.20 }; // % even row0

  // pixelav: 4.9% xtlk
  /*
  Float_t rowav[n] = { 1.654, 1.533, 1.445, 1.369, 1.325, 1.292, 1.278, 1.261 }; // rows
  Float_t eveav[n] = { 66.75, 61.41, 57.69, 54.84, 53.42, 52.49, 52.20, 51.89 }; // % even row0
  */
  // pixelav: 5.1% xtlk
  /*
  Float_t rowav[n] = { 1.691, 1.561, 1.463, 1.379, 1.332, 1.296, 1.282, 1.265 }; // rows
  Float_t eveav[n] = { 68.47, 62.75, 58.60, 55.31, 53.71, 52.66, 52.35, 51.96 }; // % even row0
  */
  // pixelav: 5.6% xtlk

  Float_t rowav[n] = { 1.787, 1.642, 1.521, 1.411, 1.353, 1.309, 1.295, 1.272 }; // rows
  Float_t eveav[n] = { 73.37, 66.84, 61.41, 56.88, 54.70, 53.31, 52.89, 52.38 }; // % even row0

  // pixelav: 6.3% xtlk fit at 1500: too steep
  /*
  Float_t rowav[n] = { 1.910, 1.764, 1.617, 1.468, 1.388, 1.332, 1.311, 1.286 }; // rows
  Float_t eveav[n] = { 79.43, 72.95, 66.30, 59.73, 56.47, 54.47, 53.77, 53.07 }; // % even row0
  */
  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // row vs thr:

  hr = new 
    TH1F( "hr",
	  "cluster size;threshold [e];cluster size [rows]",
	  16, 0, 1550 ); // axis range
  hr->SetStats(kFALSE); // no statistics
  hr->SetMinimum(1);
  hr->SetMaximum(2);
  hr->Draw();

  pr = new TGraph( n, thr, rowav );
  pr->SetLineColor(6);
  pr->SetMarkerColor(6);
  pr->SetMarkerStyle(21);
  pr->SetMarkerSize(1.5);
  pr->Draw("Pc"); // without axis option: overlay

  gr = new TGraph( n, thr, row );
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(21);
  gr->SetMarkerSize(1.5);
  gr->Draw("P"); // without axis option: overlay

  c1.Print( "even.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // even vs thr:

  he = new 
    TH1F( "he",
	  "even rowmin;threshold [e];even row cluster starts [%]",
	  16, 0, 1550 ); // axis range
  he->SetStats(kFALSE); // no statistics
  he->SetMinimum(0);
  he->SetMaximum(100);
  he->Draw();

  pe = new TGraph( n, thr, eveav );
  pe->SetLineColor(6);
  pe->SetMarkerColor(6);
  pe->SetMarkerStyle(22);
  pe->SetMarkerSize(1.5);
  pe->Draw("Pc"); // without axis option: overlay

  ge = new TGraph( n, thr, eve );
  ge->SetMarkerColor(4);
  ge->SetMarkerStyle(22);
  ge->SetMarkerSize(1.5);
  ge->Draw("P"); // without axis option: overlay

  c1.Print( "even.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "even.ps]" ); // ] closes file
  int ierr;
  ierr = system("ps2pdf even.ps");
  ierr = system("rm -f  even.ps");
  cout << "evince even.pdf" << endl;

}
