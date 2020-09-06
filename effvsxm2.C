
// Daniel Pitzl, Feb 2019
// 512 Lin eff vs x mod 200
// root -l effvsxm2.C

{
  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.022, "y" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 1.7, "y" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align
  gStyle->SetTitleX( 0.15 ); // global title
  gStyle->SetTitleY( 0.99 ); // global title

  gStyle->SetLineWidth(1);// frames
  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  gStyle->SetHistMinimumZero(); // no zero suppression

  bool ieee = 1;
  if( ieee )
    gStyle->SetOptDate(0);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //              topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 0, 0, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "effvsxm2.pdf[", "pdf" ); // [ opens file

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd;
  if( ieee )
    lgnd = new TLegend( 0.75, 0.4, 0.90, 0.8 );
  else
    lgnd = new TLegend( 0.73, 0.3, 0.93, 0.8 );

  // turn 33:

  TFile * f33 = TFile::Open( "scopeRD35580.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetTitle( "RD53A (512) Lin, 5#upoint10^{15 }p/cm^{2}, 770^{ }V" );
  effvsxm2->GetXaxis()->SetTitle( "track position within 2 pixels [#mum]" );
  effvsxm2->SetLineColor(1);
  effvsxm2->SetLineWidth(3);
  effvsxm2->SetMinimum(0.6);
  effvsxm2->SetMaximum(1.0);
  effvsxm2->SetNdivisions( -504, "y" );
  effvsxm2->Draw();
  if( ieee )
    lgnd->AddEntry( effvsxm2, " 33#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn 33#circ", "l" );

  // turn 27

  TFile * f27 = TFile::Open( "scopeRD35595.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(51); // violet
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, " 27#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn 27#circ", "l" );

  // turn 22

  TFile * f22 = TFile::Open( "scopeRD35594.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(15); // grey
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, " 22#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn 22#circ", "l" );

  // turn 17

  TFile * f17 = TFile::Open( "scopeRD35593.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(92); // gold
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, " 17#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn 17#circ", "l" );

  // turn 12

  TFile * f12 = TFile::Open( "scopeRD35592.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(8);
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, " 12#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn 12#circ", "l" );

  // turn 7

  TFile * f07 = TFile::Open( "scopeRD35591.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(4);
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, "   7#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn  7#circ", "l" );

  // turn 0

  TFile * f02 = TFile::Open( "scopeRD35596.root" );
  effvsxm2->SetStats(0);
  effvsxm2->SetLineColor(2);
  effvsxm2->SetLineWidth(3);
  effvsxm2->Draw( "same" );
  if( ieee )
    lgnd->AddEntry( effvsxm2, "   0#circ", "l" );
  else
    lgnd->AddEntry( effvsxm2, " turn  0#circ", "l" );

  lgnd->Draw( "same" );

  c1.Print( "effvsxm2.pdf" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "effvsxm2.pdf]" ); // ] closes file
  cout << "evince effvsxm2.pdf" << endl;

}
