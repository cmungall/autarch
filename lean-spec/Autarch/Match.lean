import Autarch.IR

namespace Autarch

abbrev Env := List (String × Key)

def lookup (env : Env) (name : String) : Option Key :=
  match env.find? (fun entry => entry.fst = name) with
  | some entry => some entry.snd
  | none => none

def atomMatches (env : Env) (participant : Participant) (atom : Atom) :
    Option Env :=
  match atom with
  | .exact key =>
      if participant.key = key then
        some env
      else
        none
  | .var name =>
      match lookup env name with
      | some key =>
          if participant.key = key then
            some env
          else
            none
      | none =>
          some ((name, participant.key) :: env)

def locationMatches (participant : Participant) (location : Option String) :
    Bool :=
  match location with
  | some expected => participant.location = some expected
  | none => true

def termMatches (env : Env) (participant : Participant) (term : PatTerm) :
    Option Env :=
  match atomMatches env participant term.atom with
  | some env' =>
      if locationMatches participant term.location then
        some env'
      else
        none
  | none => none

def isOptional (term : PatTerm) : Bool :=
  term.count.lo = 0

def consumeFirstMatch (env : Env) (term : PatTerm) :
    List Participant -> Option (Env × List Participant)
  | [] => none
  | participant :: participants =>
      match termMatches env participant term with
      | some env' => some (env', participants)
      | none =>
          match consumeFirstMatch env term participants with
          | some (env', remaining) =>
              some (env', participant :: remaining)
          | none => none

def matchTerms (env : Env) (participants : List Participant) :
    List PatTerm -> Option (Env × List Participant)
  | [] => some (env, participants)
  | term :: terms =>
      if isOptional term then
        match consumeFirstMatch env term participants with
        | some (env', remaining) => matchTerms env' remaining terms
        | none => matchTerms env participants terms
      else
        match consumeFirstMatch env term participants with
        | some (env', remaining) => matchTerms env' remaining terms
        | none => none

def looseMatches (reaction : Reaction) (pattern : PatReaction) : Bool :=
  match matchTerms [] reaction.lhs pattern.lhs with
  | some (env, _) =>
      match matchTerms env reaction.rhs pattern.rhs with
      | some _ => true
      | none => false
  | none => false

def strictMatches (reaction : Reaction) (pattern : PatReaction) : Bool :=
  match matchTerms [] reaction.lhs pattern.lhs with
  | some (env, lhsRemaining) =>
      match matchTerms env reaction.rhs pattern.rhs with
      | some (_, rhsRemaining) =>
          lhsRemaining.isEmpty && rhsRemaining.isEmpty
      | none => false
  | none => false

def reactionMatches (strict : Bool) (reaction : Reaction)
    (pattern : PatReaction) : Bool :=
  if strict then strictMatches reaction pattern else looseMatches reaction pattern

theorem strictImpliesLoose (reaction : Reaction) (pattern : PatReaction) :
    strictMatches reaction pattern = true ->
      looseMatches reaction pattern = true := by
  unfold strictMatches looseMatches
  cases hLeft : matchTerms [] reaction.lhs pattern.lhs with
  | none =>
      simp [hLeft]
  | some leftResult =>
      cases hRight : matchTerms leftResult.fst reaction.rhs pattern.rhs with
      | none =>
          simp [hLeft, hRight]
      | some rightResult =>
          simp [hLeft, hRight]

end Autarch
