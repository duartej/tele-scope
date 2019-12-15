
// Daniel Pitzl, Jun 2018
// RD53A resolution vs threshold
// root -l resvsthr.C

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

  gStyle->SetTitleOffset( 1.5, "x" );
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
  //                topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );

  c1.Print( "resvsthr.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  gPad->Update();// required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // RD53A 504 Lin:

  const Int_t n = 7;

  // run            33348 33470 33471 33472 33480 33481
  Float_t cur[n] = {   10,   50,   40,   30,   20,   60 }; // KRUM_CURR_LIN
  Float_t res[n] = { 9.13, 6.94, 6.86, 6.74, 6.74, 6.98 }; // [um] dxc

  Float_t sig[n];
  for( int i = 0; i < n; ++i ) {
    sig[i] = sqrt( res[i]*res[i] - 4.7*4.7 ); // subtract track
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // col vs cur:

  hcol = new 
    TH1F( "hcol",
	  "cluster size;cureshold [DAC];cluster size [columns]",
	  98, 300, 398 ); // axis range
  hcol->SetStats(kFALSE); // no statistics
  hcol->SetMinimum(1);
  hcol->SetMaximum(2);
  hcol->Draw();

  gc = new TGraph( n, cur, col );
  gc->SetMarkerColor(4);
  gc->SetMarkerStyle(22);
  gc->SetMarkerSize(1.5);
  gc->Draw("P"); // without axis option: overlay

  TLegend * lgc = new TLegend( 0.3, 0.2, 0.9, 0.3 );
  lgc->AddEntry( gc, "504 RD53A 50x150, turn 17^{o}", "p" );
  lgc->Draw( "same" );

  c1.Print( "resvscur.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // res vs cur:

  hcur = new 
    TH1F( "hcur",
	  "pixel resolution;cureshold [DAC];hit resolution [#mum]",
	  98, 300, 398 ); // axis range
  hcur->SetStats(kFALSE); // no statistics
  hcur->SetMinimum( 0);
  hcur->SetMaximum(10);
  hcur->Draw();

  gr = new TGraph( n, cur, sig );
  gr->SetMarkerColor(4);
  gr->SetMarkerStyle(20);
  gr->SetMarkerSize(1.5);
  gr->Draw("P"); // without axis option: overlay

  TLegend * lgnd = new TLegend( 0.3, 0.2, 0.9, 0.3 );
  lgnd->AddEntry( gr, "504 RD53A 50x150, turn 17^{o}", "p" );
  lgnd->Draw( "same" );

  c1.Print( "resvscur.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "resvscur.ps]" ); // ] closes file
  int ierr;
  ierr = system("ps2pdf resvscur.ps");
  ierr = system("rm -f  resvscur.ps");
  cout << "acroread resvscur.pdf" << endl;

}
