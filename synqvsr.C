
// Daniel Pitzl, Dec 2018
// 509 Syn Q0 vs R
// root -l synqvsr.C

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

  //gStyle->SetOptDate();
 
  gROOT->ForceStyle();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // square canvas:
  //               topleft x, y, width x, y
  TCanvas c1( "c1", "c1", 635, 246, 813, 837 );
  //                 to get fCw 800 fCh 800

  c1.Print( "synqvsr.ps[", "Portrait" ); // [ opens file

  gStyle->SetPaperSize( 18, 27 );

  c1.SetBottomMargin(0.15);
  c1.SetLeftMargin(0.15);
  c1.SetRightMargin(0.05);

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

  TLegend * lgnd = new TLegend( 0.83, 0.2, 0.97, 0.80 );

  // 800 V:

  TFile * f800 = TFile::Open( "scopeRD35150.root" );
  synqvsr->SetStats(0);
  synqvsr->SetTitle( "RD53 Syn 8#upoint10^{15 }p/cm^{2} " );
  synqvsr->GetXaxis()->SetTitle( "distance from irradiation center [mm]" );
  synqvsr->SetLineColor(1);
  synqvsr->SetMinimum(0);
  synqvsr->SetMaximum(18);
  synqvsr->Draw();
  lgnd->AddEntry( synqvsr, " 800 V", "l" );

  // 600 V:

  TFile * f600 = TFile::Open( "scopeRD35155.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(2);
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 600 V", "l" );

  // 500 V:

  TFile * f500 = TFile::Open( "scopeRD35156.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(4);
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 500 V", "l" );

  // 400 V:

  TFile * f400 = TFile::Open( "scopeRD35157.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(8);
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 400 V", "l" );

  // 300 V:

  TFile * f300 = TFile::Open( "scopeRD35160.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(6);
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 300 V", "l" );

  // 250 V:

  TFile * f250 = TFile::Open( "scopeRD35161.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(7); // green
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 250 V", "l" );

  // 200 V:

  TFile * f200 = TFile::Open( "scopeRD35162.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(94); // orange
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 200 V", "l" );

  // 150 V:

  TFile * f150 = TFile::Open( "scopeRD35163.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(28); // braun
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 150 V", "l" );

  // 100 V:

  TFile * f100 = TFile::Open( "scopeRD35164.root" );
  synqvsr->SetStats(0);
  synqvsr->SetLineColor(9); // blue
  synqvsr->Draw( "same" );
  lgnd->AddEntry( synqvsr, " 100 V", "l" );

  lgnd->Draw( "same" );

  c1.Print( "synqvsr.ps" );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done:

  c1.Print( "synqvsr.ps]" ); // ] closes file
  int ierr;
  ierr = system( "ps2pdf synqvsr.ps" );
  ierr = system( "rm -f  synqvsr.ps" );
  cout << "evince synqvsr.pdf" << endl;

}
