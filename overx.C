
// Feb 2018
// scopes.cc: overlay trix for alignment check
// rlq scopes34500.root
// .x overx.C
{
  gStyle->SetOptStat(111);
  //gStyle->SetStatY(0.90);
  //gStyle->SetStatW(0.30);
  //gStyle->SetStatH(0.30);

  trixc->Draw();
  N0 = trixc->Integral();
  gPad->SetGridx();

  trixclk->Draw("esame");
  N1 = trixclk->Integral();
  cout << double(N1)/N0 << endl;
}
