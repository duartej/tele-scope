
// Daniel Pitzl, May 2019
// telescope efficiency plots

// root -l tele35862.root
// .x teleff.C

{
  // quer canvas:
  //               topleft x, y, width x, y
  TCanvas * c1 = new TCanvas( "c1", "c1", 0, 0, 1600, 900 );

  // set styles:

  gStyle->SetTextFont(62);//62 = Helvetica bold
  gStyle->SetTextAlign(11);

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.02, "y");
  gStyle->SetTickLength( -0.01, "z");

  gStyle->SetLabelOffset( 0.022, "xyz" );
  gStyle->SetLabelFont( 62, "XYZ" );

  gStyle->SetTitleOffset( 1.3, "x" );
  gStyle->SetTitleOffset( 2.2, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );
  gStyle->SetTitleFont( 62, "XYZ" );

  gStyle->SetTitleBorderSize(0); // no frame around global title
  gStyle->SetTitleAlign(13); // 13 = left top align

  gStyle->SetTitleY( 0.60 ); // global title

  gStyle->SetLineWidth(2); // frames

  gStyle->SetHistLineColor(4); // 4=blau
  gStyle->SetHistLineWidth(3);
  gStyle->SetHistFillColor(5); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle(1001); // 1001 = solid

  gStyle->SetFrameLineWidth(2);

  // statistics box:

  gStyle->SetOptStat(110);
  gStyle->SetStatFormat("6.4g"); // more digits, default is 6.4g
  gStyle->SetStatFont(42); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize(1); // no 'shadow'

  gStyle->SetStatX(0.98);
  gStyle->SetStatY(0.50);
  gStyle->SetStatW(0.40);
  gStyle->SetStatH(0.35);

  gStyle->SetPalette(1); // rainbow colors

  //gStyle->SetHistMinimumZero(); // no zero suppression

  //gStyle->SetOptDate();

  c1->SetTopMargin(0.08);
  c1->SetBottomMargin(0.13);
  c1->SetLeftMargin(0.13);
  c1->SetRightMargin(0.03);

  c1->Print( "teleff.ps[", "Landscape" ); // [ opens file

  gStyle->SetPaperSize( 27, 18 );

  gROOT->ForceStyle();

  c1->Divide( 3, 2, 0, 0 ); // zero margin

  gPad->Update();

  c1->cd(1); // top left
  //gPad->SetTopMargin(0.12);
  //gPad->SetLeftMargin(0.18);
  gStyle->SetTitleX( 0.62 ); // global title
  gStyle->SetTitleSize(0.24);
  gPad->SetTickx(); // no ticks
  gPad->SetTicks(0,0); // no opposite ticks
  eff1vsx->SetLabelOffset( 0.03, "y" );
  eff1vsx->SetLabelSize( 0.05, "y" );
  eff1vsx->GetYaxis()->SetTitle("efficiency");
  eff1vsx->GetYaxis()->SetTitleSize(0.05);
  eff1vsx->SetLineWidth(3);
  eff1vsx->SetMinimum(0.8);
  eff1vsx->SetMaximum(1.0);
  eff1vsx->Draw();
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->cd(2); // top mid
  //gPad->SetTopMargin(0.12);
  gStyle->SetTitleX( 0.52 ); // global title
  gPad->SetTickx();
  gPad->SetTicky();
  gPad->SetTicks(0,0); // no opposite ticks
  eff2vsx->SetLineWidth(3);
  eff2vsx->SetMinimum(0.8);
  eff2vsx->SetMaximum(1.0);
  eff2vsx->Draw();
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->cd(3); // top right
  //gPad->SetTopMargin(0.12);
  gStyle->SetTitleX( 0.52 ); // global title
  //gPad->SetRightMargin(0.03);
  //gPad->SetTickx();
  //gPad->SetTicky();
  gPad->SetTicks(0,0); // no opposite ticks
  eff3vsx->SetLineWidth(3);
  eff3vsx->SetMinimum(0.8);
  eff3vsx->SetMaximum(1.0);
  eff3vsx->Draw();
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->cd(4); // bottom left
  //gPad->SetLeftMargin(0.18);
  gStyle->SetTitleX( 0.62 ); // global title
  //gPad->SetBottomMargin(0.15);
  gPad->SetTicks(0,0); // no opposite ticks
  eff4vsx->SetLabelOffset( 0.03, "y" );
  eff4vsx->SetLabelSize( 0.05, "y" );
  eff4vsx->GetYaxis()->SetTitle("efficiency");
  eff4vsx->GetYaxis()->SetTitleSize(0.05);
  eff4vsx->SetLabelSize( 0.05, "x" );
  eff4vsx->GetXaxis()->SetTitleSize(0.05);
  eff4vsx->GetXaxis()->CenterTitle(1);
  eff4vsx->SetLineWidth(3);
  eff4vsx->SetMinimum(0.8);
  eff4vsx->SetMaximum(1.0);
  eff4vsx->Draw(); // with axes
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->cd(5); // bottom mid
  gStyle->SetTitleX( 0.52 ); // global title
  //gPad->SetBottomMargin(0.15);
  gPad->SetTicky();
  gPad->SetTicks(0,0); // no opposite ticks
  eff5vsx->SetLabelSize( 0.05, "x" );
  eff5vsx->GetXaxis()->SetTitleSize(0.05);
  eff5vsx->GetXaxis()->CenterTitle(1);
  eff5vsx->SetLineWidth(3);
  eff5vsx->SetMinimum(0.8);
  eff5vsx->SetMaximum(1.0);
  eff5vsx->Draw(); // ticksx
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->cd(6); // bottom right
  gStyle->SetTitleX( 0.52 ); // global title
  //gPad->SetRightMargin(0.03);
  //gPad->SetBottomMargin(0.15);
  gPad->SetTicky();
  gPad->SetTicks(0,0); // no opposite ticks
  eff6vsx->SetLabelSize( 0.05, "x" );
  eff6vsx->GetXaxis()->SetTitleSize(0.05);
  eff6vsx->GetXaxis()->CenterTitle(1);
  eff6vsx->SetLineWidth(3);
  eff6vsx->SetMinimum(0.8);
  eff6vsx->SetMaximum(1.0);
  eff6vsx->Draw(); // ticksx
  cout << gPad->GetWNDC() << "  " << gPad->GetHNDC() << endl;

  c1->Print( "teleff.ps" );

  c1->Print( "teleff.ps]" ); // ] closes file

  int retval = 0;
  retval = system( "ps2pdf teleff.ps" );
  retval = system( "rm -f  teleff.ps" );
  cout << "evince teleff.pdf" << endl;

}
