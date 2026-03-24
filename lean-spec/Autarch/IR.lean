namespace Autarch

inductive Key where
  | chebi : String -> Key
  | smiles : String -> Key
  | name : String -> Key
  | anonymous : String -> Key
deriving Repr, DecidableEq, Inhabited

structure Participant where
  key : Key
  count : Nat := 1
  location : Option String := none
  stoichiometry : Option String := none
  polymerIndex : Option String := none
deriving Repr, DecidableEq, Inhabited

structure Reaction where
  lhs : List Participant
  rhs : List Participant
  label : String := ""
  transport : Bool := false
deriving Repr, DecidableEq, Inhabited

inductive Atom where
  | exact : Key -> Atom
  | var : String -> Atom
deriving Repr, DecidableEq, Inhabited

structure Count where
  lo : Nat := 1
  hi : Option Nat := some 1
deriving Repr, DecidableEq, Inhabited

structure PatTerm where
  atom : Atom
  count : Count := {}
  location : Option String := none
deriving Repr, DecidableEq, Inhabited

structure PatReaction where
  lhs : List PatTerm
  rhs : List PatTerm
  label : String := ""
deriving Repr, DecidableEq, Inhabited

end Autarch
