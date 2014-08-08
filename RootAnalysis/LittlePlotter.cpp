#include "LittlePlotter.h"
//TH1F wrangling code:
void LittlePlotter::book(TH1F* hBase)
{
	std::vector<TH1F*> vHists;
	for(std::vector<TString>::const_iterator catIt = m_categories.begin(); catIt != m_categories.end(); ++catIt)
	{
		TString cat = *catIt;
		TH1F* hCat = (TH1F*) hBase->Clone(cat+TString("_")+hBase->GetName());
		hCat->SetTitle(cat);//Store category as Title for ease of access
		vHists.push_back(hCat);
	}
	m_th1f.insert(std::pair<TString, std::vector<TH1F*> >(hBase->GetName(), vHists));
	delete hBase;
}
TH1F* LittlePlotter::getPlot(TString& plotName, TString& category)
{
	TH1F* selectedPlot;
	mIt_th1f = m_th1f.find(plotName);	
	if(mIt_th1f != m_th1f.end())
	{
		for(std::vector<TH1F*>::iterator histIt = mIt_th1f->second.begin(); histIt != mIt_th1f->second.end(); ++histIt)
		{
			TH1F* hist = *histIt;
			if(hist->GetTitle() == category)
			{
				return hist;	
			}
		}
		std::cout<<"LittlePlotter::getPlot: Could not find histogram "<< plotName <<" for "<< category << std::endl;
		printAllPlotNames();
		return NULL;
	}else{
		std::cout<<"LittlePlotter::getPlot: Could not find histogram "<< plotName <<" for "<< category << std::endl;
		printAllPlotNames();
		return NULL;
	}	
}
std::vector<TH1F*> LittlePlotter::getPlots(TString& plotName, std::vector<TString>& categories)
{
	std::vector<TH1F*> selectedPlots;
	mIt_th1f = m_th1f.find(plotName);	
	if(mIt_th1f != m_th1f.end())
	{
		for(std::vector<TH1F*>::iterator histIt = mIt_th1f->second.begin(); histIt != mIt_th1f->second.end(); ++histIt)
		{
			TH1F* hist = *histIt;
			for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
			{
				if(hist->GetTitle() == *cat)
				{
					selectedPlots.push_back(hist);
				}
			}
		}
	}else{
		std::cout<<"Could not find histogram "<< plotName << std::endl;
	}	
	return selectedPlots;
}
void LittlePlotter::scalePlots(TString plotName, std::vector<TString>& categories ,std::map<TString,float>& numerator, std::map<TString,float>& denominator)
{
	mIt_th1f = m_th1f.find(plotName);	
	if(mIt_th1f != m_th1f.end())
	{
		for(std::vector<TH1F*>::iterator histIt = mIt_th1f->second.begin(); histIt != mIt_th1f->second.end(); ++histIt)
		{
			TH1F* hist = *histIt;
			for(std::vector<TString>::iterator cat = categories.begin(); cat != categories.end(); ++cat)
			{
				if(hist->GetTitle() == *cat)
				{
				  std::map<TString, float>::const_iterator numIt = numerator.find(*cat);
				  float num = numIt->second;
				  std::map<TString, float>::const_iterator denomIt = denominator.find(*cat);
				  float denom = denomIt->second;
				  hist->Scale(num/denom);
				}
			}
		}
	}else{
		std::cout<<"Could not find histogram "<< plotName << std::endl;
	}	
	return;
}
//Plotting code:
void LittlePlotter::plotAlone(TString plotName, std::vector<TString>& categories)
{
	std::vector<TH1F*> plots = getPlots(plotName, categories);
	if(plots.size() < 1) return;
	TCanvas* c0 = new TCanvas("c0", "Overlay plots", 800, 600);
	for(size_t i = 0; i < plots.size(); ++i)
	{
		TH1F* plot = plots[i];
		//float effSig = effSigma(plot);
		//std::cout<<"LittlePlotter: "<< plot->GetName() <<", effSigma = "<< effSig << std::endl;
		TPaveText* txt = new TPaveText(0.65, 0.5, 0.9, 0.9, "NDC");
		txt->SetBorderSize(0);
		txt->SetFillColor(0);
		//txt->AddText(Form("Effective #sigma = %.3f GeV", effSig));
		//txt->AddText(Form("RMS = %.3f GeV", plot->GetRMS()));
		//txt->AddText(Form("Effective #sigma/<m_{b#bar{b}}> = %.3f", effSig/plot->GetMean()));
		//txt->AddText(Form("RMS/<m_{b#bar{b}}> = %.3f", plot->GetRMS()/plot->GetMean()));
		plot->SetLineWidth(2);
		plot->Draw("E1");			
		txt->Draw();
		c0->SaveAs(TString("Plots/")+plot->GetName()+TString(".pdf"));
		TH1F* plotNorm = (TH1F*) plot->Clone(TString(plot->GetName())+"Norm");
		plotNorm->Scale(1./plot->Integral());
		plotNorm->Draw("HIST");			
		txt->Draw();
		plotNorm->GetYaxis()->SetTitle("Arbitrary Units");
		c0->SaveAs(TString("Plots/")+plot->GetName()+TString("_norm.pdf"));
		delete txt;
	}	
	delete c0;
} 
void LittlePlotter::plotOverlay(TString plotName, std::vector<TString>& categories)
{
	std::vector<TH1F*> plots = getPlots(plotName, categories);
	if(plots.size() < 1) return;
	TCanvas* cOver = new TCanvas("cOver", "Overlay plots", 800, 600);
	TCanvas* cOverNorm = new TCanvas("cOverNorm", "Normalised Overlay plots", 800, 600);
	TLegend leg(0.72, 0.72, 0.92, 0.92);
	//TLegend leg(0.60, 0.60, 0.90, 0.90);
	leg.SetBorderSize(0);
	leg.SetFillColor(0);
	float normMax = 0.; TH1F* firstNorm;
	for(size_t i = 0; i < plots.size(); ++i)
	{
		TH1F* plot = plots[i];
		plot->SetLineColor(m_sigCol[i]);
		plot->SetLineStyle(i+1);
		plot->SetLineWidth(2);
		TH1F* plotNorm = (TH1F*) plot->Clone(TString(plot->GetName())+"Norm");
		plotNorm->Scale(1./plot->Integral());
		if(plotNorm->GetMaximum() > normMax) normMax = plotNorm->GetMaximum();
		if(i == 0)
		{
			firstNorm = plotNorm;
			cOver->cd();
			plot->Draw("HIST");			
			cOverNorm->cd();
			plotNorm->Draw("HIST");			
		}else{
			cOver->cd();
			plot->Draw("HISTsame");	
			cOverNorm->cd();
			plotNorm->Draw("HISTsame");	
		}
		leg.AddEntry(plot, legendName(plot->GetTitle()), "L");
	}	
	cOver->cd();
	leg.Draw();
	cOver->SaveAs(TString("Plots/")+plots[0]->GetName()+TString(".pdf"));
	cOverNorm->cd();
	firstNorm->GetYaxis()->SetTitle("Arbitrary Units");
	firstNorm->GetYaxis()->SetRangeUser(0, 1.05*normMax);
	leg.Draw();
	cOverNorm->SaveAs(TString("Plots/")+plots[0]->GetName()+TString("_norm.pdf"));
	//cOverNorm->SaveAs(TString("Plots/")+plots[0]->GetName()+TString("_norm.eps"));
	delete cOver; delete cOverNorm;
} 
void LittlePlotter::plotSoBVsEff(TString plotName, std::vector<TString>& signals, std::vector<TString>& backgrounds)
{
	std::vector<TH1F*> sigPlots = getPlots(plotName, signals);
	std::vector<TH1F*> backPlots = getPlots(plotName, backgrounds);

	TCanvas* cSoB = new TCanvas("cSoB", "S/B Plots", 800, 600);
	TLegend leg(0.60, 0.70, 0.90, 0.90);
	leg.SetBorderSize(0);
	leg.SetFillColor(0);

	TH1F* sumSig = sigPlots[0]; 
	sumSig->SetLineColor(m_sigCol[0]);
	sumSig->SetLineStyle(1);
	sumSig->SetLineWidth(3);
	for(size_t i = 1; i < sigPlots.size(); ++i)
	{
		TH1F* plot = sigPlots[i];
		sumSig->Add(plot);
	}	
	TH1F* sumBack = backPlots[0];
	sumBack->SetLineColor(kBlack);
	sumBack->SetFillColor(m_backCol[0]);
	for(size_t i = 1; i < backPlots.size(); ++i)
	{
		TH1F* plot = backPlots[i];
		sumBack->Add(plot);
	}	
	int nBinsData = 0;
	for(int i = 0; i < sumSig->GetNbinsX()+1; ++i)
	{
		double sig =sumSig->Integral(i, sumSig->GetNbinsX()+1); 
		double back =sumBack->Integral(i, sumBack->GetNbinsX()+1); 
		if(sig > 1e-4 && back > 1e-4) ++nBinsData;
	}
	
	//TVectorD eff(sumSig->GetNbinsX()+1), soB(sumSig->GetNbinsX()+1), soSqrtB(sumSig->GetNbinsX()+1);
	TVectorD eff(nBinsData), soB(nBinsData), soSqrtB(nBinsData);
	double nSig = sumSig->Integral(0, sumSig->GetNbinsX()+1);
	//for(int i = 0; i < sumSig->GetNbinsX()+1; ++i)
	for(int i = 0; i < nBinsData; ++i)
	{
		double sig =sumSig->Integral(i, sumSig->GetNbinsX()+1); 
		double back =sumBack->Integral(i, sumBack->GetNbinsX()+1); 
		eff[i] = sig/nSig;
		soB[i] = sig/back;
		soSqrtB[i] = sig/sqrt(back);
		//std::cout<<"eff = "<< eff[i] <<", S/B = "<< soB[i] <<", S/√B = "<< soSqrtB[i] << std::endl;
	}
	TGraph gSoB(eff, soB);
	gSoB.SetMarkerStyle(21);
	gSoB.Draw("AP");
	gSoB.GetXaxis()->SetTitle("Signal Efficiency");
	gSoB.GetYaxis()->SetTitle("S/B");
	cSoB->SaveAs(TString("Plots/")+sigPlots[0]->GetName()+TString("_SoB.pdf"));
	TGraph gSoSqrtB(eff, soSqrtB);
	gSoSqrtB.SetMarkerStyle(21);
	gSoSqrtB.Draw("AP");
	gSoSqrtB.GetXaxis()->SetTitle("Signal Efficiency");
	gSoSqrtB.GetYaxis()->SetTitle("S/#sqrt{B}");
	cSoB->SaveAs(TString("Plots/")+sigPlots[0]->GetName()+TString("_SoSqrtB.pdf"));
	delete cSoB;
}
void LittlePlotter::plotStack(TString plotName, std::vector<TString>& signals, std::vector<TString>& backgrounds, bool sigAtBottom)
{
	std::vector<TH1F*> sigPlots = getPlots(plotName, signals);
	std::vector<TH1F*> backPlots = getPlots(plotName, backgrounds);

	TCanvas* cStack = new TCanvas("cStack", "Stacked Signal & Background Plots", 800, 600);
	TLegend leg(0.60, 0.70, 0.90, 0.90);
	leg.SetBorderSize(0);
	leg.SetFillColor(0);

	for(size_t i = 0; i < sigPlots.size(); ++i)
	{
		TH1F* plot = sigPlots[i];
		plot->SetLineColor(m_sigCol[i]);
		plot->SetLineStyle(i+1);
		plot->SetLineWidth(3);
		leg.AddEntry(plot, legendName(plot->GetTitle()), "L");
	}	
	THStack backStack("backStack", "Background");
	for(size_t i = 0; i < backPlots.size(); ++i)
	{
		TH1F* plot = backPlots[i];
		plot->SetLineColor(kBlack);
		plot->SetFillColor(m_backCol[i]);
		leg.AddEntry(plot, legendName(plot->GetTitle()), "F");
		backStack.Add(plot);
	}	
	//std::vector<THStack*> backStacks;
	//if(sigAtBottom) {
	//  if(backStack.GetMaximum() > sigPlots[0]->GetMaximum()) ;
	//  backStack.Draw("HIST");
	//  backStack.GetXaxis()->SetTitle(sigPlots[0]->GetXaxis()->GetTitle());
	//  backStack.GetYaxis()->SetTitle(sigPlots[0]->GetYaxis()->GetTitle());	
	//  for(size_t i = 0; i < sigPlots.size(); ++i)
	//    {
	//      sigPlots[i]->Draw("sameHIST");
	//    }	
	//}
	if(sigAtBottom) {
	  if(backStack.GetMaximum() > sigPlots[0]->GetMaximum()) {
	    backStack.Draw("HIST");
	    backStack.GetXaxis()->SetTitle(sigPlots[0]->GetXaxis()->GetTitle());
	    backStack.GetYaxis()->SetTitle(sigPlots[0]->GetYaxis()->GetTitle());	
	  }
	  else {
	    sigPlots[0]->Draw("HIST");
	    sigPlots[0]->GetXaxis()->SetTitle(sigPlots[0]->GetXaxis()->GetTitle());
	    sigPlots[0]->GetYaxis()->SetTitle(sigPlots[0]->GetYaxis()->GetTitle());	
	  }
	  backStack.Draw("sameHIST");
	  for(size_t i = 0; i < sigPlots.size(); ++i)
	    {
	      sigPlots[i]->Draw("sameHIST");
	    }	
	}
	else {
	  for(size_t i = 0; i < sigPlots.size(); ++i)
	    {
	      TH1F* plot = sigPlots[i];
	      THStack* newStack = (THStack*) backStack.Clone();
	      newStack->Add(plot);
		if(i==0) {
		  newStack->Draw("HIST");
		  newStack->GetXaxis()->SetTitle(sigPlots[0]->GetXaxis()->GetTitle());
		  newStack->GetYaxis()->SetTitle(sigPlots[0]->GetYaxis()->GetTitle());
		}
		else newStack->Draw("sameHIST");
	    }
	  backStack.Draw("sameHIST");	
	}
	
	leg.Draw();
	TPaveText* txt3 = new TPaveText(0.65, 0.5, 0.9, 0.7, "NDC");
	txt3->SetBorderSize(0);
	txt3->SetFillColor(0);
	txt3->AddText("Expected for #intLdt = 3000 fb^{-1}");
	txt3->Draw();
	backStack.Draw("sameAXIS");
	cStack->SaveAs(TString("Plots/")+sigPlots[0]->GetName()+TString("_Stack.pdf"));
	//cStack->SaveAs(TString("Plots/")+sigPlots[0]->GetName()+TString("_Stack.eps"));
	delete cStack;
}
std::pair<float, float> LittlePlotter::calculateSignificance(TString plotName, TString signal, std::vector<TString>& backgrounds, float minVal, float maxVal)
{
        std::cout << "" << std::endl;
	std::cout << "LittlePlotter:calculateSignificance for " << signal << " in mass range " << minVal << " - " << maxVal << std::endl;
	TH1F* sigPlot = getPlot(plotName, signal);
	int minBin = sigPlot->FindBin(minVal);
	int maxBin = sigPlot->FindBin(maxVal);
	double sError = 0.; double bError = 0.;
	float s = sigPlot->IntegralAndError(minBin, maxBin, sError);
	std::cout<<"LittlePlotter:calculateSignificance: S = "<< s << " +/- " << sError << std::endl;
	std::vector<TH1F*> backPlots = getPlots(plotName, backgrounds);
	TH1F* totBack = new TH1F("totBack", "", backPlots[0]->GetNbinsX(), backPlots[0]->GetXaxis()->GetXmin(), backPlots[0]->GetXaxis()->GetXmax());
	for(size_t i = 0; i < backPlots.size(); ++i)
	{
		totBack->Add(backPlots[i]);	
	}
	float b = totBack->IntegralAndError(minBin, maxBin, bError);
	std::cout<<"LittlePlotter:calculateSignificance: B = "<< b << " +/- " << bError << std::endl;
	float sqrtB = sqrt(b);
	float significance = s/sqrtB;

	//float error = sqrt((1/b)*(bError*bError + (s*s*sError*sError)/(4*b*b)));
	float error = (sqrt(((sError*sError)/(s*s)) + ((bError*bError)/(b*b*4.))))*significance;
	std::cout<<"LittlePlotter::calculateSignificance: S/sqrt(B) = "<< significance << " +/-" << error << std::endl;
	delete totBack;
	return std::pair<float, float>(significance, error);
}
void LittlePlotter::setupColourPalettes()
{
	gStyle->SetOptTitle(0);
	gStyle->SetOptStat(0);
	gStyle->SetOptFit(0);
	// put tick marks on top and RHS of plots
	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
	m_sigCol[0] = kGreen+2; m_sigCol[3] = kMagenta+4; m_sigCol[2] = kBlue; m_sigCol[1] = kBlue-6;
	m_sigCol[4] = kAzure+4; m_sigCol[5] = kCyan+3; m_sigCol[6] = kTeal; m_sigCol[7] = kTeal-7;
	m_backCol[0] = kRed+1; m_backCol[1] = kRed-8; m_backCol[2] = kOrange+2; m_backCol[3] = kOrange+7;
	m_backCol[4] = kOrange-1; m_backCol[5] = kOrange-2; m_backCol[6] = kOrange-3; m_backCol[7] = kOrange-4;
}

TString LittlePlotter::legendName(const char* plotTitle)
{
  TString pLab = TString(plotTitle);
  if(pLab.BeginsWith("RSG")) { 
    pLab.Remove(0,3); pLab.Insert(0,"G_{KK} (M = "); pLab.Append(" GeV)"); 
    return pLab;
  }
  if(pLab.EqualTo("TTbar")) {
    pLab = "t#bar{t}";
    return pLab;
  } 
  if(pLab.EqualTo("Sherpa")) {
    pLab = "QCD Multi-jet";
    return pLab;
  } 
  return pLab;
}
void LittlePlotter::writeAll()
{
	for(mIt_th1f = m_th1f.begin(); mIt_th1f != m_th1f.end(); ++mIt_th1f)
	{
		for(std::vector<TH1F*>::iterator hIt = mIt_th1f->second.begin(); hIt != mIt_th1f->second.end(); ++hIt)
		{
			TH1F* h = *hIt;
			h->Write();
		}
	}
}
void LittlePlotter::printAllPlotNames()
{
	std::cout<<"LittlePlotter: existing plots are: "<<std::endl;
	for(mIt_th1f = m_th1f.begin(); mIt_th1f != m_th1f.end(); ++mIt_th1f)
	{
		for(std::vector<TH1F*>::iterator hIt = mIt_th1f->second.begin(); hIt != mIt_th1f->second.end(); ++hIt)
		{
			TH1F* h = *hIt;
			std::cout<< h->GetName() << std::endl;
		}
	}

}
void LittlePlotter::printAllCategories()
{
	std::cout<<"LittlePlotter: defined event categories are: "<< std::endl;
	for(std::vector<TString>::iterator cat = m_categories.begin(); cat != m_categories.end(); ++cat)
	{
		std::cout<< *cat << std::endl;
	}
}
