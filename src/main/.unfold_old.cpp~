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
  mode_ptr mode = organism.getMode();
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
  
  ptree& af_node = paramlist_node.add("Param", "");
  af_node.put("<xmlattr>.name","BindingAffinity");
  af_node.put("<xmlattr>.value"     ,"");
  af_node.put("<xmlattr>.limit_low" ,1);
  af_node.put("<xmlattr>.limit_high",1);
  af_node.put("<xmlattr>.tweak"     ,0);
  
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
      quench_node.put("<xmlattr>.value",      -coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  -coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", -coef_ptrs[0]->getLimHigh());
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
      quench_node.put("<xmlattr>.value",      -coefs[0]);
      quench_node.put("<xmlattr>.limit_low",  -coef_ptrs[0]->getLimLow());
      quench_node.put("<xmlattr>.limit_high", -coef_ptrs[0]->getLimHigh());
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
      drep_node.put("<xmlattr>.limit_low",  0);
      drep_node.put("<xmlattr>.limit_high", 0);
      drep_node.put("<xmlattr>.tweak",      0);
      
      ptree& dd_node = lparamlist_node.add("Param","");
      dd_node.put("<xmlattr>.name", "DirectRepDistance");
      dd_node.put("<xmlattr>.value", 100);
      dd_node.put("<xmlattr>.limit_low",  50);
      dd_node.put("<xmlattr>.limit_high", 500);
      dd_node.put("<xmlattr>.tweak",      0);
      
      ptree& ddd_node = lparamlist_node.add("Param","");
      ddd_node.put("<xmlattr>.name",       "DirectRepDistanceDelta");
      ddd_node.put("<xmlattr>.value",      50);
      ddd_node.put("<xmlattr>.limit_low",  50);
      ddd_node.put("<xmlattr>.limit_high", 50);
      ddd_node.put("<xmlattr>.tweak",      0);
      
      ptree& act_node = lparamlist_node.add("Param","");
      act_node.put("<xmlattr>.name",       "ActivationCoef");
      act_node.put("<xmlattr>.value",      coefs[1]);
      act_node.put("<xmlattr>.limit_low",  coef_ptrs[1]->getLimLow());
      act_node.put("<xmlattr>.limit_high", coef_ptrs[1]->getLimHigh());
      act_node.put("<xmlattr>.tweak",      (int) coef_ptrs[1]->isAnnealed());
    
    }
    
    ptree& kmax_node = lparamlist_node.add("Param","");
    double_param_ptr kmax_param = tf.getKmaxParam();
    kmax_node.put("<xmlattr>.name",       "k_max");
    kmax_node.put("<xmlattr>.value",      kmax_param->getValue());
    kmax_node.put("<xmlattr>.limit_low",  kmax_param->getLimLow());
    kmax_node.put("<xmlattr>.limit_high", kmax_param->getLimHigh());
    kmax_node.put("<xmlattr>.tweak",      (int) kmax_param->isAnnealed());
    
    ptree& lambda_node = lparamlist_node.add("Param","");
    double_param_ptr lambda_param = tf.getLambdaParam();
    lambda_node.put("<xmlattr>.name",       "lambda");
    lambda_node.put("<xmlattr>.value",      lambda_param->getValue());
    lambda_node.put("<xmlattr>.limit_low",  lambda_param->getLimLow());
    lambda_node.put("<xmlattr>.limit_high", lambda_param->getLimHigh());
    lambda_node.put("<xmlattr>.tweak",      (int) lambda_param->isAnnealed());
    
    ptree& threshold_node = lparamlist_node.add("Param","");
    double_param_ptr threshold_param = tf.getThresholdParam();
    threshold_node.put("<xmlattr>.name",       "threshold");
    threshold_node.put("<xmlattr>.value",      threshold_param->getValue());
    threshold_node.put("<xmlattr>.limit_low",  threshold_param->getLimLow());
    threshold_node.put("<xmlattr>.limit_high", threshold_param->getLimHigh());
    threshold_node.put("<xmlattr>.tweak",      (int) threshold_param->isAnnealed());
    
    ptree& weightmatrix_node = ligand_node.add("WeightMatrix","");
    weightmatrix_node.put("<xmlattr>.name", tf.getName());
    weightmatrix_node.put("<xmlattr>.bsize", tf.getBindingSize());
    weightmatrix_node.put("<xmlattr>.pseudo", 1);
    weightmatrix_node.put("<xmlattr>.type", "Align");
    
    pwm_param_ptr pwmparam = tf.getPWMParam();
    weightmatrix_node.put("<xmlattr>.maxscore", pwmparam->getValue().getMaxScore());
    weightmatrix_node.put("<xmlattr>.minscore", 0);
    vector<vector<double> > pwm = pwmparam->getValue().getPWM(PCM);
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
    double_param_ptr kcoopparam = coops[i]->getKcoopParam();
    kcoop.put("<xmlattr>.name",       "Kcoop");
    kcoop.put("<xmlattr>.value",      kcoopparam->getValue());
    kcoop.put("<xmlattr>.tweak",      (int) kcoopparam->isAnnealed());
    kcoop.put("<xmlattr>.limit_high", kcoopparam->getLimHigh());
    kcoop.put("<xmlattr>.limit_low",  kcoopparam->getLimLow());
    
    map<string, double_param_ptr>& distparams = coops[i]->getDist()->getParams();
    double_param_ptr distA = distparams["A"];
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name",       "CoopDistance");
    cd.put("<xmlattr>.value",      distA->getValue());
    cd.put("<xmlattr>.tweak",      (int) distA->isAnnealed());
    cd.put("<xmlattr>.limit_high", distA->getLimHigh());
    cd.put("<xmlattr>.limit_low",  distA->getLimLow());
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoopDistanceDelta");
    cdd.put("<xmlattr>.value", 10);
    cdd.put("<xmlattr>.tweak", 0);
    cdd.put("<xmlattr>.limit_high", 10);
    cdd.put("<xmlattr>.limit_low", 10);
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
    
    double_param_ptr eff = coeffects[i]->getEfficiencyParam();
    ptree& cparamlist = coop_node.add("ParamList","");
    ptree& kcoop = cparamlist.add("Param","");
    kcoop.put("<xmlattr>.name", "CoactCoef");
    kcoop.put("<xmlattr>.value", coeffects[i]->getEfficiency());
    kcoop.put("<xmlattr>.tweak", (int) eff->isAnnealed());
    kcoop.put("<xmlattr>.limit_high", eff->getLimHigh());
    kcoop.put("<xmlattr>.limit_low", eff->getLimLow());
    
    map<string, double_param_ptr>& distparams = coeffects[i]->getDist()->getParams();
    double_param_ptr distA = distparams["A"];
    double_param_ptr distB = distparams["B"];
    
    ptree& cd = cparamlist.add("Param","");
    cd.put("<xmlattr>.name", "CoactDistance");
    cd.put("<xmlattr>.value",      distA->getValue());
    cd.put("<xmlattr>.tweak",      (int) distA->isAnnealed());
    cd.put("<xmlattr>.limit_high", distA->getLimHigh());
    cd.put("<xmlattr>.limit_low",  distA->getLimLow());
    
    ptree& cdd = cparamlist.add("Param","");
    cdd.put("<xmlattr>.name", "CoactDistanceDelta");
    cdd.put("<xmlattr>.value",      distB->getValue());
    cdd.put("<xmlattr>.tweak",      (int) distB->isAnnealed());
    cdd.put("<xmlattr>.limit_high", distB->getLimHigh());
    cdd.put("<xmlattr>.limit_low",  distB->getLimLow());
  }
  
  ptree& mode_node = paramset_node.add("Mode","");
  mode_node.add("<xmlattr>.dynamic",1);
  mode_node.add("<xmlattr>.overlappingQuench",0);
  mode_node.add("<xmlattr>.overlappingCoact",0);
  mode_node.add("<xmlattr>.directRepression",0);
  mode_node.add("<xmlattr>.cooperativity",1);
  mode_node.add("<xmlattr>.scaledChisq",1);
  mode_node.add("<xmlattr>.showGroups",0);
  mode_node.add("<xmlattr>.appcomp",0);
  mode_node.add("<xmlattr>.SmoothTxnrate",1);
  mode_node.add("<xmlattr>.hilde",0);
  mode_node.add("<xmlattr>.random",0);
  
  table_ptr ratedata = organism.getRateData();
  ptree& constructlist_node = paramset_node.add("ConstructList","");
  int ngenes = genes->size();
  vector<string>& lins = ratedata->getRowNames();
  for (int i = 0; i<ngenes; i++)
  {
    ptree& construct_node = constructlist_node.add("Construct","");
    Gene& gene = genes->getGene(i);
    construct_node.put("<xmlattr>.genotype", gene.getName());
    construct_node.put("<xmlattr>.include",  (int) gene.getInclude());
    
    scale_factor_ptr scale = gene.getScale();
    double_param_ptr scaleA = scale->getA();
    ptree& scale_node = construct_node.add("PositionEffect","");
    scale_node.put("<xmlattr>.scale_factor",scaleA->getValue());
    scale_node.put("<xmlattr>.limit_low",   scaleA->getLimLow());
    scale_node.put("<xmlattr>.limit_high",  scaleA->getLimHigh());
    scale_node.put("<xmlattr>.tweak",       (int) scaleA->isAnnealed());
    
    ptree& seq_node = construct_node.add("Sequence","");
    seq_node.put("<xmlattr>.left_bound", gene.getLeftBound());
    seq_node.put("<xmlattr>.right_bound", gene.getRightBound());
    seq_node.put("<xmlattr>.AT_bkgd", 1-mode->getGC()/2);
    seq_node.put("<xmlattr>.CG_bkgd", mode->getGC()/2);
    seq_node.put("<xmlattr>.tweak", (int) gene.getSequenceParam()->isAnnealed());
    seq_node.put("", gene.getSequenceString());
    
    ptree& sitelist_node = construct_node.add("BindingSiteList","");
    
    vector<double*>& col = ratedata->getCol(gene.getName());
    ptree& conclist_node = construct_node.add("MrnaConcList","");
    conclist_node.put("<xmlattr>.timepoint", 85.175000);
    int nrow = lins.size();
    for (int j=0; j<nrow; j++)
    {
      ptree& conc_node = conclist_node.add("MrnaConc","");
      conc_node.put("<xmlattr>.lin",  lins[j]);
      conc_node.put("<xmlattr>.conc", *(col[j]));
    }
      
    
    
  }
  ptree& data_node = system_node.add("Data","");
  ptree& ligandconclist_node = data_node.add("LigandConcList","");
  ligandconclist_node.put("<xmlattr>.timepoint", 85.175000);
  int nlins = lins.size();
  table_ptr tfdata = organism.getTFData();
  vector<string>& tfnames = tfdata->getColNames();
  int ntfdata = tfnames.size();
  for (int i=0; i<nlins; i++)
  {
    ptree& ligandconc_node = ligandconclist_node.add("LigandConc","");
    ligandconc_node.put("<xmlattr>.lin", lins[i]);
    for (int j=0; j<ntfdata; j++)
    {
      string name = tfnames[j];
      string id = m[name];
      ligandconc_node.put("<xmlattr>."+id, tfdata->getDataPoint("ID",lins[i],"TF",tfnames[j]));
      //cerr << "ID: " << lins[i] << endl;
      //cerr << "TF: " << tfnames[j] << endl;
      //cerr << "data: " << tfdata->getDataPoint("ID",lins[i],"TF",tfnames[j]) << endl;
    }
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


