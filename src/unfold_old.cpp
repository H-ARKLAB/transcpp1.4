/*********************************************************************************
*                                                                                *
*     unfold.cpp                                                                 *
*                                                                                *
*     Reads output and allows user to print various data                         *
*                                                                                *
*********************************************************************************/
#include "pwm.h"
#include "TF.h"
#include "gene.h"

#include "datatable.h"
#include "twobit.h"
#include "organism.h"

#include <fstream>
#include <boost/property_tree/xml_parser.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/version.hpp>
#include <boost/lexical_cast.hpp>

#include <unistd.h>
#include <getopt.h>

using boost::property_tree::ptree;
#define to_string_ boost::lexical_cast<string>


int mode_verbose;

map<string,string> ligand_map()
{
  map<string,string> m;
  m["bcd"]   = "B";
  m["cad"]   = "C";
  m["dst"]   = "D";
  m["hb"]    = "H";
  m["Kr"]    = "K";
  m["kni"]   = "N";
  m["gt"]    = "G";
  m["tll"]   = "T";
  m["slp"]   = "S";
  m["prd"]   = "P";
  m["eve"]   = "E";
  m["run"]   = "R";
  m["odd"]   = "O";
  m["ftz"]   = "F";
  m["dic"]   = "I";
  m["hairy"] = "Y";
  m["moo"]   = "M";
  m["zld"]   = "Z";
  return m;
}
  

void print_xml(Organism& organism, string& section)
{
  ptree out;
  ptree& system_node = out.add("System","");
  
  ptree& description_node = system_node.add("Description", "converted from new transcpp code");
  
  ptree& paramset_node = system_node.add("ParamSet","");
  if (section == "Input")
    paramset_node.put("<xmlattr>.section", "input");
  else
    paramset_node.put("<xmlattr>.section", "eqparms");
  
  nuclei_ptr nuclei = organism.getNuclei();
  int nnuc = nuclei->size();
  
  ptree& paramlist_node = paramset_node.add("ParamList","");
  ptree& nnuc_node = paramlist_node.add("Param", "");
  nnuc_node.put("<xmlattr>.name","NoNuclei");
  nnuc_node.put("<xmlattr>.value",nnuc);
  
  ptree& fgf_node = paramlist_node.add("Param", "");
  fgf_node.put("<xmlattr>.name","Fgf");
  fgf_node.put("<xmlattr>.value"     ,1);
  fgf_node.put("<xmlattr>.limit_low" ,1);
  fgf_node.put("<xmlattr>.limit_high",1);
  fgf_node.put("<xmlattr>.tweak"     ,0);
  
  ptree& rmax_node = paramlist_node.add("Param", "");
  rmax_node.put("<xmlattr>.name","MaxRate");
  rmax_node.put("<xmlattr>.value",     255);
  rmax_node.put("<xmlattr>.limit_low" ,255);
  rmax_node.put("<xmlattr>.limit_high",255);
  rmax_node.put("<xmlattr>.tweak"     ,0);
  
  genes_ptr    genes    = organism.getGenes();
  promoter_ptr promoter = genes->getGene(0).getPromoter();
  map<string, double_param_ptr>& promoter_params = promoter->getParamMap();
  ptree& theta_node = paramlist_node.add("Param", "");
  theta_node.put("<xmlattr>.name",       "ActivationThreshold");
  theta_node.put("<xmlattr>.value",      promoter_params["Theta"]->getValue());
  theta_node.put("<xmlattr>.limit_low",  promoter_params["Theta"]->getLimLow());
  theta_node.put("<xmlattr>.limit_high", promoter_params["Theta"]->getLimHigh());
  theta_node.put("<xmlattr>.tweak",      (int) promoter_params["Theta"]->isAnnealed());
  
  ptree& q_node = paramlist_node.add("Param", "");
  q_node.put("<xmlattr>.name",       "Q");
  q_node.put("<xmlattr>.value",      promoter_params["Q"]->getValue());
  q_node.put("<xmlattr>.limit_low",  promoter_params["Q"]->getLimLow());
  q_node.put("<xmlattr>.limit_high", promoter_params["Q"]->getLimHigh());
  q_node.put("<xmlattr>.tweak",      (int) promoter_params["Q"]->isAnnealed());
  
  tfs_ptr tfs = organism.getTFs();
  int ntfs = tfs->size();
  map<string,string> m = ligand_map();
  
  ptree& ligandlist_node = paramset_node.add("LigandList","");
  
  for (int i=0; i<ntfs; i++)
  {
    TF& tf = tfs->getTF(i);
    vector<double> coefs = tf.getCoefs();
    vector<double_param_ptr>& coef_ptrs = tf.getCoefPtrs();
    
    int ligand_type;
    if (coefs.size() == 2)
      ligand_type = 2;
    else if (coefs[0] >= 0)
      ligand_type = 0;
    else if (coefs[0] < 0)
      ligand_type = 1;
    
    ptree& ligand_node = ligandlist_node.add("Ligand","");
  
    ligand_node.put("<xmlattr>.name", tf.getName());
    ligand_node.put("<xmlattr>.ligand_id", m[tf.getName()]);
    if (coefs.size() == 2)
      ligand_node.put("<xmlattr>.ligand_type", 2);
    else if (coefs[0] >= 0)
      ligand_node.put("<xmlattr>.ligand_type", 0);
    else if (coefs[0] < 0)
      ligand_node.put("<xmlattr>.ligand_type", 1);
    ligand_node.put("<xmlattr>.include", 1);
    
    ptree& lparamlist_node = ligand_node.add("ParamList","");
    
    if (ligand_type == 0)
    {
      ptree& act_node = lparamlist_node.add("Param","");
      act_node.put("<xmlattr>.name",       "ActivationCoef");
      act_node.put("<xmlattr>.value",      coefs[0]);
      act_node.put("<xmlattr>.limit_low",  coef_ptrs[0]->getLimLow());
      act_node.put("<xmlattr>.limit_high", coef_ptrs[0]->getLimHigh());
      act_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());

    } 
    else if (ligand_type == 1)
    {
      ptree& quench_node = lparamlist_node.add("Param","");
      quench_node.put("<xmlattr>.name", "QuenchingCoef");
      quench_node.put("<xmlattr>.value",      coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", coef_ptrs[0]->getLimHigh());
      quench_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());
      
      ptree& qd_node = lparamlist_node.add("Param","");
      qd_node.put("<xmlattr>.name",       "QuenchingDistance");
      qd_node.put("<xmlattr>.value",      100);
      qd_node.put("<xmlattr>.limit_low",  100);
      qd_node.put("<xmlattr>.limit_high", 100);
      qd_node.put("<xmlattr>.tweak",      0);
      
      ptree& qdd_node = lparamlist_node.add("Param","");
      qdd_node.put("<xmlattr>.name",       "QuenchingDistanceDelta");
      qdd_node.put("<xmlattr>.value",      50);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 50);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& drep_node = lparamlist_node.add("Param","");
      drep_node.put("<xmlattr>.name",      "DirectRepressionCoef");
      drep_node.put("<xmlattr>.value",      0);
      drep_node.put("<xmlattr>.limit_low",  0);
      drep_node.put("<xmlattr>.limit_high", 0);
      drep_node.put("<xmlattr>.tweak",      0);
      
      ptree& dd_node = lparamlist_node.add("Param","");
      dd_node.put("<xmlattr>.name",       "DirectRepDistance");
      dd_node.put("<xmlattr>.value",      100);
      dd_node.put("<xmlattr>.limit_low",  100);
      dd_node.put("<xmlattr>.limit_high", 500);
      dd_node.put("<xmlattr>.tweak",      0);
      
      ptree& ddd_node = lparamlist_node.add("Param","");
      ddd_node.put("<xmlattr>.name",       "DirectRepDistanceDelta");
      ddd_node.put("<xmlattr>.value",      50);
      ddd_node.put("<xmlattr>.limit_low",  50);
      ddd_node.put("<xmlattr>.limit_high", 50);
      ddd_node.put("<xmlattr>.tweak",      0);
    }
    else if (ligand_type == 2)
    {
      ptree& quench_node = lparamlist_node.add("Param","");
      quench_node.put("<xmlattr>.name",       "QuenchingCoef");
      quench_node.put("<xmlattr>.value",      coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", coef_ptrs[0]->getLimHigh());
      quench_node.put("<xmlattr>.tweak",      (int) coef_ptrs[0]->isAnnealed());
      
      ptree& qd_node = lparamlist_node.add("Param","");
      qd_node.put("<xmlattr>.name",       "QuenchingDistance");
      qd_node.put("<xmlattr>.value",      100);
      qd_node.put("<xmlattr>.limit_low",  100);
      qd_node.put("<xmlattr>.limit_high", 100);
      qd_node.put("<xmlattr>.tweak",      0);
      
      ptree& qdd_node = lparamlist_node.add("Param","");
      qdd_node.put("<xmlattr>.name",       "QuenchingDistanceDelta");
      qdd_node.put("<xmlattr>.value",      50);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 50);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& drep_node = lparamlist_node.add("Param","");
      drep_node.put("<xmlattr>.name",      "DirectRepressionCoef");
      drep_node.put("<xmlattr>.value",     0);
      qdd_node.put("<xmlattr>.limit_low",  0);
      qdd_node.put("<xmlattr>.limit_high", 0);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& dd_node = lparamlist_node.add("Param","");
      dd_node.put("<xmlattr>.name", "DirectRepDistance");
      dd_node.put("<xmlattr>.value", 100);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 500);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& ddd_node = lparamlist_node.add("Param","");
      ddd_node.put("<xmlattr>.name",       "DirectRepDistanceDelta");
      ddd_node.put("<xmlattr>.value",      50);
      qdd_node.put("<xmlattr>.limit_low",  50);
      qdd_node.put("<xmlattr>.limit_high", 50);
      qdd_node.put("<xmlattr>.tweak",      0);
      
      ptree& act_node = lparamlist_node.add("Param","");
      act_node.put("<xmlattr>.name",       "ActivationCoef");
      act_node.put("<xmlattr>.value",      coefs[1]);
      qdd_node.put("<xmlattr>.limit_low",  coef_ptrs[1]->getLimLow());
      qdd_node.put("<xmlattr>.limit_high", coef_ptrs[1]->getLimHigh());
      qdd_node.put("<xmlattr>.tweak",      (int) coef_ptrs[1]->isAnnealed());
    
    }
    
    ptree& kmax_node = lparamlist_node.add("Param","");
    kmax_node.put("<xmlattr>.name", "k_max");
    kmax_node.put("<xmlattr>.value", tf.getKmax());
    
    ptree& lambda_node = lparamlist_node.add("Param","");
    lambda_node.put("<xmlattr>.name", "lambda");
    lambda_node.put("<xmlattr>.value", tf.getLambda());
    
    ptree& threshold_node = lparamlist_node.add("Param","");
    threshold_node.put("<xmlattr>.name", "threshold");
    threshold_node.put("<xmlattr>.value", tf.getThreshold());
    
    ptree& weightmatrix_node = ligand_node.add("WeightMatrix","");
    weightmatrix_node.put("<xmlattr>.name", tf.getName());
    weightmatrix_node.put("<xmlattr>.bsize", tf.getBindingSize());
    weightmatrix_node.put("<xmlattr>.pseudo", 1);
    
    vector<vector<double> > pwm = tf.getPWMParam()->getValue().getPWM(PCM);
    int pwmlen = pwm.size();
    for (int j=0; j<pwmlen; j++)
    {
      ptree& position_node = weightmatrix_node.add("Position", "");
      position_node.put("<xmlattr>.Weights", to_string_(round(pwm[j][0])) + ";"
                                           + to_string_(round(pwm[j][1])) + ";"
                                           + to_string_(round(pwm[j][2])) + ";"
                                           + to_string_(round(pwm[j][3])));
    }
  }
  ptree& gintlist_node = paramset_node.add("GlobalInteractionList","");
  vector<coop_ptr>& coops = organism.getCoops()->getAllCoops();
  int ncoops = coops.size();
  for (int i=0; i<ncoops; i++)
  {
    ptree& coop_node = gintlist_node.add("Coop","");
    coop_node.put("<xmlattr>.actor", m[coops[i]->getTFs().first]);
    coop_node.put("<xmlattr>.target", m[coops[i]->getTFs().second]);
    coop_node.put("<xmlattr>.mode", "pair-wise");
    coop_node.put("<xmlattr>.action_type", 1);
    coop_node.put("<xmlattr>.include", 1);
    
    ptree& cparamlist = coop_node.add("ParamList","");
    ptree& kcoop = cparamlist.add("Param","");
    kcoop.put("<xmlattr>.name", "Kcoop");
    kcoop.put("<xmlattr>.value", coops[i]->getK());
    
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name", "CoopDistance");
    cd.put("<xmlattr>.value", coops[i]->getDist()->getParams()["Max"]);
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoopDistanceDelta");
    cdd.put("<xmlattr>.value", 10);
  }
  
  vector<coeffect_ptr>& coeffects = organism.getCoeffects()->getAllCoeffects();
  int ncoeffects = coeffects.size();
  for (int i=0; i<ncoeffects; i++)
  {
    ptree& coop_node = gintlist_node.add("Coact","");
    coop_node.put("<xmlattr>.actor", m[coeffects[i]->getActor()]);
    coop_node.put("<xmlattr>.target", m[coeffects[i]->getTarget()]);
    coop_node.put("<xmlattr>.action_type", 0);
    coop_node.put("<xmlattr>.include", 1);
    
    ptree& cparamlist = coop_node.add("ParamList","");
    ptree& kcoop = cparamlist.add("Param","");
    kcoop.put("<xmlattr>.name", "CoactCoef");
    kcoop.put("<xmlattr>.value", coeffects[i]->getEfficiency());
    
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name", "CoactDistance");
    cd.put("<xmlattr>.value", coeffects[i]->getDist()->getParams()["A"]);
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoactDistanceDelta");
    cdd.put("<xmlattr>.value", coeffects[i]->getDist()->getParams()["A"]);
  }
  
#if BOOST_VERSION / 100 % 1000 < 56
  write_xml_element(cout, 
    basic_string<ptree::key_type::value_type>(), 
    out, 
    -1, 
    boost::property_tree::xml_writer_make_settings<char>(' ', 2));
#else
  write_xml_element(cout, 
    basic_string<ptree::key_type::value_type>(), 
    out, 
    -1, 
    boost::property_tree::xml_writer_make_settings<string>(' ', 2));
#endif
  
}
    
    
void print_f(Organism& organism, string& gname)
{
  ptree out;
  ptree& occ_node = out.add("Occupancy","");
  occ_node.put("<xmlattr>.genotype", gname);
  
}
  
static const char *optString = "fhosXx:";

static const struct option longOpts[] = {
    { "help",        no_argument,       NULL, 'h' },
    { "section",     required_argument, NULL, 'x' },
    { "occupancy",   no_argument,       NULL, 'f' },
    { "sequence",    no_argument,       NULL, 's' },
    { "outputxml",   no_argument,       NULL, 'o' }
    //{ 0, 0, 0, 0}
};

void display_usage()
{
  cerr << endl << "\t Usage" << endl << endl
       << "\t unfold [options] input_file genotype" << endl << endl
       << "\t Options" << endl
       << "\t --help      [-h]  print this message" << endl
       << "\t --section   [-x]  use section of input file (default eqparms)" << endl
       << "\t --sequence  [-s]  use sequence-level code" << endl
       << "\t --occupancy [-f]  prints binding site fractional occupancy" << endl
       << "\t --outputxml [-o]  outputs the old xml format" << endl << endl;
  exit(1);
}

int main(int argc, char* argv[])
{
  int opt = 0;
  int longIndex = 0;
  string section_name("eqparms");
  
  bool sequence  = false;
  bool occupancy = false;
  bool outputxml = false;
  bool xmlflag   = false;

  string infile_name;
  
  opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  while(opt != -1)
  {
    switch (opt)
    {
      case 'f':
        occupancy = true;
        break;
      case 'h':
        display_usage();
        break;
      case 'o':
        outputxml = true;
        break;
      case 's':
        sequence = true;
        break;
      case 'X':
        xmlflag = true;
        break;
      case 'x':
        section_name = optarg;
        break;
      case 0: 
        display_usage();
        break;
      case '?':
        cerr << "?" << endl;
        display_usage();
        break;
    }
    opt = getopt_long( argc, argv, optString, longOpts, &longIndex );
  }
  
  // convert to new node names
  if (section_name == "eqparms")
    section_name = "Output";
  else if (section_name == "input")
    section_name = "Input";
  
  if (argc <= optind)
    error("Missing input file");
  
  infile_name = argv[optind];
  
  if (infile_name == "")
    display_usage();
  
  ifstream infile(infile_name.c_str());

  ptree pt;
  read_xml(infile, pt, boost::property_tree::xml_parser::trim_whitespace);
  
  ptree& root_node   = pt.get_child("Root");
  ptree& mode_node   = root_node.get_child("Mode");
  ptree& section_node = root_node.get_child(section_name);
  
  mode_ptr mode(new Mode(infile_name,mode_node));

  mode->setVerbose(0);
  Organism embryo(section_node, mode);
  
  if (outputxml)
    print_xml(embryo,section_name);
  
  return 0;
}


