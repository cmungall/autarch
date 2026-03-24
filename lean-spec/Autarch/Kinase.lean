import Autarch.Match

namespace Autarch

def ATP : Key := .chebi "CHEBI:30616"
def ADP : Key := .chebi "CHEBI:456216"
def GTP : Key := .chebi "CHEBI:37565"
def GDP : Key := .chebi "CHEBI:58189"
def Proton : Key := .chebi "CHEBI:15378"

def optionalProton : PatTerm := {
  atom := .exact Proton,
  count := { lo := 0, hi := some 1 }
}

def kinaseATP : PatReaction := {
  lhs := [
    { atom := .exact ATP },
    { atom := .var "substrate" },
    optionalProton
  ],
  rhs := [
    { atom := .exact ADP },
    { atom := .var "product" },
    optionalProton
  ],
  label := "ATP kinase pattern"
}

def kinaseGTP : PatReaction := {
  lhs := [
    { atom := .exact GTP },
    { atom := .var "substrate" },
    optionalProton
  ],
  rhs := [
    { atom := .exact GDP },
    { atom := .var "product" },
    optionalProton
  ],
  label := "GTP kinase pattern"
}

def kinasePatterns : List PatReaction := [kinaseATP, kinaseGTP]

def Glucose : Key := .chebi "CHEBI:17234"
def Glucose6Phosphate : Key := .chebi "CHEBI:17665"

def simpleKinaseReaction : Reaction := {
  lhs := [
    { key := ATP },
    { key := Glucose }
  ],
  rhs := [
    { key := ADP },
    { key := Glucose6Phosphate }
  ],
  label := "ATP + glucose -> ADP + glucose-6-phosphate"
}

example : looseMatches simpleKinaseReaction kinaseATP = true := by
  decide

example : strictMatches simpleKinaseReaction kinaseATP = true := by
  decide

theorem kinaseStrictImpliesLoose :
    strictMatches simpleKinaseReaction kinaseATP = true ->
      looseMatches simpleKinaseReaction kinaseATP = true := by
  exact strictImpliesLoose simpleKinaseReaction kinaseATP

end Autarch
