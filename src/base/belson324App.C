#include "belson324App.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
belson324App::validParams()
{
  InputParameters params = MooseApp::validParams();
  params.set<bool>("use_legacy_material_output") = false;
  params.set<bool>("use_legacy_initial_residual_evaluation_behavior") = false;
  return params;
}

belson324App::belson324App(InputParameters parameters) : MooseApp(parameters)
{
  belson324App::registerAll(_factory, _action_factory, _syntax);
}

belson324App::~belson324App() {}

void
belson324App::registerAll(Factory & f, ActionFactory & af, Syntax & syntax)
{
  ModulesApp::registerAllObjects<belson324App>(f, af, syntax);
  Registry::registerObjectsTo(f, {"belson324App"});
  Registry::registerActionsTo(af, {"belson324App"});

  /* register custom execute flags, action syntax, etc. here */
}

void
belson324App::registerApps()
{
  registerApp(belson324App);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
belson324App__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  belson324App::registerAll(f, af, s);
}
extern "C" void
belson324App__registerApps()
{
  belson324App::registerApps();
}
