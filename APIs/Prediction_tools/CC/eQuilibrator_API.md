**eQuilibrator API Reference**

# equilibrator_api.model.bounds

Define lower and upper bounds on compounds.

**Module Contents**
```python
class equilibrator_api.model.bounds.BaseBounds
    Bases: object
    A base class for declaring bounds on things.

    abstract copy()
        Return a (deep) copy of self.

    abstract get_lower_bound(compound: Union[str, equilibrator_cache.Compound])
        Get the lower bound for this key.
            Parameters:
                key – a compound

    abstract get_upper_bound(compound: Union[str, equilibrator_cache.Compound])
        Get the upper bound for this key.
            Parameters:
                key – a compound

    get_lower_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Iterable[equilibrator_api.Q_]
        Get the bounds for a set of keys in order.
        Parameters:

            compounds – an iterable of Compounds or strings
        Returns:
            an iterable of the lower bounds

    get_upper_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Iterable[equilibrator_api.Q_]
        Get the bounds for a set of keys in order.
        Parameters:
            compounds – an iterable of Compounds or strings
        Returns:
            an iterable of the upper bounds

    get_bound_tuple(compound: Union[str, equilibrator_cache.Compound]) → Tuple[equilibrator_api.Q_, equilibrator_api.Q_]
        Get both upper and lower bounds for this key.
        Parameters:
            compound – a Compound object or string
        Returns:
            a 2-tuple (lower bound, upper bound)

    get_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Tuple[Iterable[equilibrator_api.Q_], Iterable[equilibrator_api.Q_]]
        Get the bounds for a set of compounds.
        Parameters:
            compounds – an iterable of Compounds
        Returns:
            a 2-tuple (lower bounds, upper bounds)

    static conc2ln_conc(b: equilibrator_api.Q_) → float
        Convert a concentration to log-concentration.
        Parameters:
            b – a concentration
        Returns:
            the log concentration

    get_ln_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Tuple[Iterable[float], Iterable[float]]
        Get the log-bounds for a set of compounds.
        Parameters:
            compounds – an iterable of Compounds or strings
        Returns:
            a 2-tuple (log lower bounds, log upper bounds)

    get_ln_lower_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Iterable[float]
        Get the log lower bounds for a set of compounds.
        Parameters:
            compounds – an iterable of Compounds or strings
        Returns:
            an iterable of log lower bounds

    get_ln_upper_bounds(compounds: Iterable[Union[str, equilibrator_cache.Compound]]) → Iterable[float]
        Get the log upper bounds for a set of compounds.
        Parameters:
            compounds – an iterable of Compounds or strings
        Returns:
            an iterable of log upper bounds

    set_bounds(compound: Union[str, equilibrator_cache.Compound], lb: equilibrator_api.Q_, ub: equilibrator_api.Q_) → None
        Set bounds for a specific key.
        Parameters:
            key – a Compounds or string
            lb – the lower bound value
            ub – the upper bound value
```
```python
class equilibrator_api.model.bounds.Bounds(lower_bounds: Dict[Union[str, equilibrator_cache.Compound], equilibrator_api.Q_] = None, upper_bounds: Dict[Union[str, equilibrator_cache.Compound], equilibrator_api.Q_] = None, default_lb: equilibrator_api.Q_ = default_conc_lb, default_ub: equilibrator_api.Q_ = default_conc_ub)
    Bases: BaseBounds
    Contains upper and lower bounds for various keys.
    Allows for defaults.

    DEFAULT_BOUNDS

    classmethod from_csv(f: TextIO, comp_contrib: equilibrator_api.ComponentContribution, default_lb: equilibrator_api.Q_ = default_conc_lb, default_ub: equilibrator_api.Q_ = default_conc_ub) → Bounds
        Read Bounds from a CSV file.
        Parameters:
            f (TextIO) – an open .csv file stream
            comp_contrib (ComponentContribution) – used for parsing compound accessions
            default_lb (Q_) – the default lower bound
            default_ub (Q_) – the default upper bound

    to_data_frame() → pandas.DataFrame
        Convert the list of bounds to a Pandas DataFrame.

    check_bounds() → None
        Assert the bounds are valid (i.e. that lb <= ub).

    copy() → Bounds
        Return a deep copy of self.

    get_lower_bound(compound: Union[str, equilibrator_cache.Compound]) → equilibrator_api.Q_
        Get the lower bound for this compound.

    get_upper_bound(compound: Union[str, equilibrator_cache.Compound]) → equilibrator_api.Q_
        Get the upper bound for this compound.

    static get_default_bounds(comp_contrib: equilibrator_api.ComponentContribution) → Bounds
        Return the default lower and upper bounds for a pre-determined list.
        Parameters:
            comp_contrib (ComponentContribution) –
        Returns:
            a Bounds object with the default values
```

# equilibrator_api.model.model

A basic stoichiometric model with thermodynamics.

**Module Contents**
```python
class equilibrator_api.model.model.StoichiometricModel(S: pandas.DataFrame, compound_dict: Dict[str, equilibrator_cache.Compound], reaction_dict: Dict[str, equilibrator_cache.Reaction], comp_contrib: Optional[equilibrator_api.ComponentContribution] = None, standard_dg_primes: Optional[equilibrator_api.Q_] = None, dg_sigma: Optional[equilibrator_api.Q_] = None, bounds: Optional[equilibrator_api.model.Bounds] = None, config_dict: Optional[Dict[str, str]] = None)
    Bases: object
    A basic stoichiometric model with thermodynamics.
    Designed as a base model for 'Pathway' which also includes flux directions and magnitudes.

    MINIMAL_STDEV = 0.001

    configure() → None
        Configure the Component Contribution aqueous conditions.

    property compound_ids → Iterable[str]
        Get the list of compound IDs.

    property compounds → Iterable[equilibrator_cache.Compound]
        Get the list of Compound objects.

    property compound_df → pandas.DataFrame
        Get a DataFrame with all the compound data.
        The columns are: compound_id, lower_bound, upper_bound

    property reaction_ids → Iterable[str]
        Get the list of reaction IDs.

    property reactions → Iterable[equilibrator_cache.Reaction]
        Get the list of Reaction objects.

    property reaction_formulas → Iterable[str]
        Iterate through all the reaction formulas.
        Returns:
            the reaction formulas

    property reaction_df → pandas.DataFrame
        Get a DataFrame with all the reaction data.
        The columns are: reaction_id, reaction_formula, standard_dg_prime

    update_standard_dgs() → None
        Calculate the standard G' values and uncertainties.
        Use the Component Contribution method.

    set_bounds(cid: str, lb: Optional[equilibrator_api.Q_] = None, ub: Optional[equilibrator_api.Q_] = None) → None
        Set the lower and upper bound of a compound.
        Parameters:
            compound_id (str) – the compound ID
            lb (Quantity, optional) – the new concentration lower bound (ignored if the value is None)
            ub (Quantity, optional) – the new concentration upper bound (ignored if the value is None)

    get_bounds(cid: str) → Tuple[equilibrator_api.Q_, equilibrator_api.Q_]
        Get the lower and upper bound of a compound.
        Parameters:
            compound_id (str) – the compound ID
        Returns:
            lb (Quantity, optional) – the new concentration lower bound (ignored if the value is None)
            ub (Quantity, optional) – the new concentration upper bound (ignored if the value is None)

    property bounds → Tuple[Iterable[equilibrator_api.Q_], Iterable[equilibrator_api.Q_]]
        Get the concentration bounds.
        The order of compounds is according to the stoichiometric matrix index.
        Returns:
            tuple of (lower bounds, upper bounds)

    property bound_df → pandas.DataFrame
        Get a DataFrame with all the bounds data.

    property ln_conc_lb → numpy.array
        Get the log lower bounds on the concentrations.
        The order of compounds is according to the stoichiometric matrix index.
        Returns:
            a NumPy array of the log lower bounds

    property ln_conc_ub → numpy.ndarray
        Get the log upper bounds on the concentrations.
        The order of compounds is according to the stoichiometric matrix index.
        Returns:
            a NumPy array of the log upper bounds

    property ln_conc_mu → numpy.array
        Get mean of log concentration distribution based on the bounds.
        The order of compounds is according to the stoichiometric matrix index.
        Returns:
            a NumPy array with the mean of the log concentrations

    property ln_conc_sigma → numpy.array
        Get stdev of log concentration distribution based on the bounds.
        Returns:
            a NumPy array with the stdev of the log concentrations

    static read_thermodynamics(thermo_sbtab: equilibrator_api.model.SBtabTable, config_dict: Dict[str, str]) → Dict[str, equilibrator_api.Q_]
        Read the 'thermodynamics' table from an SBtab.
        Parameters:
            thermo_sbtab (SBtabTable) – A SBtabTable containing the thermodynamic data
            config_dict (dict) – A dictionary containing the configuration arguments
        Returns:
            A dictionary mapping reaction IDs to standard G' values.

    classmethod from_network_sbtab(filename: Union[str, equilibrator_api.model.SBtabDocument], comp_contrib: Optional[equilibrator_api.ComponentContribution] = None, freetext: bool = True, bounds: Optional[equilibrator_api.model.Bounds] = None) → StoichiometricModel
        Initialize a Pathway object using a 'network'-only SBtab.
        Parameters:
            filename (str, SBtabDocument) – a filename containing an SBtabDocument (or the SBtabDocument object itself) defining the network (topology) only
            comp_contrib (ComponentContribution, optional) – a ComponentContribution object needed for parsing and searching the reactions. also used to set the aqueous parameters (pH, I, etc.)
            freetext (bool, optional) – a flag indicating whether the reactions are given as freetext (i.e. common names for compounds) or by standard database accessions (Default value: True)
            bounds (Bounds, optional) – bounds on metabolite concentrations (by default uses the "data/cofactors.csv" file in equilibrator-api)
        Returns:
            a Pathway object

    classmethod from_sbtab(filename: Union[str, equilibrator_api.model.SBtabDocument], comp_contrib: Optional[equilibrator_api.ComponentContribution] = None) → StoichiometricModel
        Parse and SBtabDocument and return a StoichiometricModel.
        Parameters:
            filename (str or SBtabDocument) – a filename containing an SBtabDocument (or the SBtabDocument object itself) defining the pathway
            comp_contrib (ComponentContribution, optional) – a ComponentContribution object needed for parsing and searching the reactions. also used to set the aqueous parameters (pH, I, etc.)
        Returns:
            stoich_model – A StoichiometricModel object based on the configuration SBtab

    to_sbtab() → equilibrator_api.model.SBtabDocument
        Export the model to an SBtabDocument.

    write_sbtab(filename: str) → None
        Write the pathway to an SBtab file.
```

```python
equilibrator_api.model.open_sbtabdoc(filename: Union[str, sbtab.SBtab.SBtabDocument]) → sbtab.SBtab.SBtabDocument
    Open a file as an SBtabDocument.
    Checks whether it is already an SBtabDocument object, otherwise reads the CSV file and returns the parsed object.
```
# equilibrator_api.compatibility

Provide functions for compatibility with COBRA.

**Module Contents**
```python
equilibrator_api.compatibility.map_cobra_reactions(cache: equilibrator_cache.CompoundCache, reactions: List[cobra.Reaction], **kwargs) → Dict[str, equilibrator_api.phased_reaction.PhasedReaction]
    Translate COBRA reactions to eQuilibrator phased reactions.
    Parameters:
        cache (equilibrator_cache.CompoundCache) –
        reactions (iterable of cobra.Reaction) – A list of reactions to map to equilibrator phased reactions.
        kwargs – Any further keyword arguments are passed to the underlying function for mapping metabolites.
    Returns:
        A mapping from COBRA reaction identifiers to equilibrator phased reactions where such a mapping can be established.
    See also:
        equilibrator_cache.compatibility.map_cobra_metabolites
```

# equilibrator_api.component_contribution

A wrapper for the GibbeEnergyPredictor in component-contribution.

**Module Contents**
```python
equilibrator_api.component_contribution.find_most_abundant_ms(cpd: equilibrator_cache.Compound, p_h: equilibrator_api.Q_, p_mg: equilibrator_api.Q_, ionic_strength: equilibrator_api.Q_, temperature: equilibrator_api.Q_) → equilibrator_cache.CompoundMicrospecies
    Find the most abundant microspecies based on transformed energies.
```
```python
equilibrator_api.component_contribution.predict_protons_and_charge(rxn: equilibrator_api.phased_reaction.PhasedReaction, p_h: equilibrator_api.Q_, p_mg: equilibrator_api.Q_, ionic_strength: equilibrator_api.Q_, temperature: equilibrator_api.Q_) → Tuple[float, float, float]
    Find the #protons and charge of a transport half-reaction.
```
```python
class equilibrator_api.component_contribution.ComponentContribution(rmse_inf: equilibrator_api.Q_ = default_rmse_inf, ccache: Optional[equilibrator_cache.CompoundCache] = None, predictor: Optional[component_contribution.predict.GibbsE] = None)
    Bases: object
    A wrapper class for GibbsEnergyPredictor.
    Also holds default conditions for compounds in the different phases.

    property p_h → equilibrator_api.Q_
        Get the pH.

    property p_mg → equilibrator_api.Q_
        Get the pMg.

    property ionic_strength → equilibrator_api.Q_
        Get the ionic strength.

    property temperature → equilibrator_api.Q_
        Get the temperature.

    static legacy() → ComponentContribution
        Initialize a ComponentContribution object with the legacy version.
        The legacy version is intended for compatibility with older versions of equilibrator api (0.2.x - 0.3.1). Starting from 0.3.2, there is a significant change in the predictions caused by an improved Mg2+ concentration model.
        Returns:
            A ComponentContribution object

    static initialize_custom_version(rmse_inf: equilibrator_api.Q_ = default_rmse_inf, ccache_settings: component_contribution.ZenodoSettings = DEFAULT_COMPOUND_CACHE_SETTINGS, cc_params_settings: component_contribution.ZenodoSettings = DEFAULT_CC_PARAMS_SETTINGS) → ComponentContribution
        Initialize ComponentContribution object with custom Zenodo versions.
        Parameters:
            rmse_inf (Quantity, optional) – Set the factor by which to multiply the error covariance matrix for reactions outside the span of Component Contribution. (Default value: 1e-5 kJ/mol)
            settings (ZenodoSettings) – The doi, filename and md5 of the
        Returns:
            A ComponentContribution object

    get_compound(compound_id: str) → Union[equilibrator_cache.Compound, None]
        Get a Compound using the DB namespace and its accession.
        Returns:
            cpd

    get_compound_by_inchi(inchi: str) → Union[equilibrator_cache.Compound, None]
        Get a Compound using InChI.
        Returns:
            cpd

    search_compound_by_inchi_key(inchi_key: str) → List[equilibrator_cache.Compound]
        Get a Compound using InChI.
        Returns:
            cpd

    search_compound(query: str) → Union[None, equilibrator_cache.Compound]
        Try to find the compound that matches the name best.
        Parameters:
            query (str) – an (approximate) compound name
        Returns:
            cpd – the best match

    parse_reaction_formula(formula: str) → equilibrator_api.phased_reaction.PhasedReaction
        Parse reaction text using exact match.
        Parameters:
            formula (str) – a string containing the reaction formula
        Returns:
            rxn

    search_reaction(formula: str) → equilibrator_api.phased_reaction.PhasedReaction
        Search a reaction written using compound names (approximately).
        Parameters:
            formula (str) – a string containing the reaction formula
        Returns:
            rxn

    balance_by_oxidation(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.phased_reaction.PhasedReaction
        Convert an unbalanced reaction into an oxidation reaction.
        By adding H2O, O2, Pi, CO2, and NH4+ to both sides.

    get_oxidation_reaction(compound: equilibrator_cache.Compound) → equilibrator_api.phased_reaction.PhasedReaction
        Generate an oxidation Reaction for a single compound.
        Generate a Reaction object which represents the oxidation reaction of this compound using O2. If there are N atoms, the product must be NH3 (and not N2) to represent biological processes. Other atoms other than C, N, H, and O will raise an exception.

    property RT → equilibrator_api.Q_
        Get the value of RT.

    standard_dg_formation(compound: equilibrator_cache.Compound) → Tuple[Optional[float], Optional[numpy.ndarray]]
        Get the (mu, sigma) predictions of a compound's formation energy.
        Parameters:
            compound (Compound) – a Compound object
        Returns:
            mu (float) – the mean of the standard formation Gibbs energy estimate
            sigma_fin (array) – a vector representing the square root of the covariance matrix (uncertainty)
            sigma_inf (array) – a vector representing the infinite-uncertainty eigenvalues of the covariance matrix

    standard_dg(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the chemical reaction energies of a reaction.
        Returns:
            standard_dg – the dG0 in kJ/mol and standard error. To calculate the 95% confidence interval, use the range -1.96 to 1.96 times this value

    standard_dg_multi(reactions: List[equilibrator_api.phased_reaction.PhasedReaction], uncertainty_representation: str = 'cov') → Tuple[numpy.ndarray, numpy.ndarray]
        Calculate the chemical reaction energies of a list of reactions.
        Using the major microspecies of each of the reactants.
        Parameters:
            reactions (List[PhasedReaction]) – a list of PhasedReaction objects to estimate
            uncertainty_representation (str) – which representation to use for the uncertainties. cov would return a full covariance matrix. precision would return the precision matrix (i.e. the inverse of the covariance matrix). sqrt would return a sqaure root of the covariance, based on the uncertainty vectors. fullrank would return a full-rank square root of the covariance which is a compressed form of the sqrt result. (Default value: cov)
        Returns:
            standard_dg (Quantity) – the estimated standard reaction Gibbs energies based on the the major microspecies
            dg_uncertainty (Quantity) – the uncertainty matrix (in either 'cov', 'sqrt' or 'fullrank' format)

    standard_dg_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the transformed reaction energies of a reaction.
        Returns:
            standard_dg – the dG0_prime in kJ/mol and standard error. To calculate the 95% confidence interval, use the range -1.96 to 1.96 times this value

    dg_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the dG'0 of a single reaction.
        Returns:
            dg – the dG_prime in kJ/mol and standard error. To calculate the 95% confidence interval, use the range -1.96 to 1.96 times this value

    standard_dg_prime_multi(reactions: List[equilibrator_api.phased_reaction.PhasedReaction], uncertainty_representation: str = 'cov', minimize_norm: bool = False) → Tuple[equilibrator_api.Q_, equilibrator_api.Q_]
        Calculate the transformed reaction energies of a list of reactions.
        Parameters:
            reactions (List[PhasedReaction]) – a list of PhasedReaction objects to estimate
            uncertainty_representation (str) – which representation to use for the uncertainties. cov would return a full covariance matrix. precision would return the precision matrix (i.e. the inverse of the covariance matrix). sqrt would return a sqaure root of the covariance, based on the uncertainty vectors. fullrank would return a full-rank square root of the covariance which is a compressed form of the sqrt result. (Default value: cov)
            minimize_norm (bool) – if True, use an orthogonal projection to minimize the norm2 of the result vector (keeping it within the finite-uncertainty sub-space, i.e. only moving along eigenvectors with infinite uncertainty).
        Returns:
            standard_dg_prime (Quantity) – the CC estimation of the reactions' standard transformed energies
            dg_uncertainty (Quantity) – the uncertainty co-variance matrix (in either 'cov', 'sqrt' or 'fullrank' format)

    physiological_dg_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the dG'm of a single reaction.
        Assume all aqueous reactants are at 1 mM, gas reactants at 1 mbar and the rest at their standard concentration.
        Returns:
            standard_dg_primes (ndarray) – a 1D NumPy array containing the CC estimates for the reactions' physiological dG'
            dg_sigma (ndarray) – the second is a 2D numpy array containing the covariance matrix of the standard errors of the estimates. one can use the eigenvectors of the matrix to define a confidence high-dimensional space, or use dg_sigma as the covariance of a Gaussian used for sampling (where 'standard_dg_primes' is the mean of that Gaussian).

    dgf_prime_sensitivity_to_p_h(compound: equilibrator_cache.Compound) → equilibrator_api.ureg.Quantity
        Calculate the sensitivity of the chemical formation energy to pH.
        Returns:
            The derivative of ∆𝐺𝑓 with respect to pH, in kJ/mol.

    dg_prime_sensitivity_to_p_h(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Quantity
        Calculate the sensitivity of the chemical reaction energy to pH.
        Returns:
            The derivative of ∆𝐺𝑟 with respect to pH, in kJ/mol.

    ln_reversibility_index(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the reversibility index (ln Gamma) of a single reaction.
        Returns:
            ln_RI – the reversibility index (in natural log scale).

    standard_e_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the E'0 of a single half-reaction.
        Returns:
            standard_e_prime (ureg.Measurement)
            the estimated standard electrostatic potential of reaction and E0_uncertainty is the standard deviation of estimation. Multiply it by 1.96 to get a 95% confidence interval (which is the value shown on eQuilibrator).

    physiological_e_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the E'0 of a single half-reaction.
        Returns:
            physiological_e_prime (ureg.Measurement)
            the estimated physiological electrostatic potential of reaction and E0_uncertainty is the standard deviation of estimation. Multiply it by 1.96 to get a 95% confidence interval (which is the value shown on eQuilibrator).

    e_prime(reaction: equilibrator_api.phased_reaction.PhasedReaction) → equilibrator_api.ureg.Measurement
        Calculate the E'0 of a single half-reaction.
        Returns:
            e_prime (ureg.Measurement)
            the estimated electrostatic potential of reaction and E0_uncertainty is the standard deviation of estimation. Multiply it by 1.96 to get a 95% confidence interval (which is the value shown on eQuilibrator).

    dg_analysis(reaction: equilibrator_api.phased_reaction.PhasedReaction) → List[Dict[str, object]]
        Get the analysis of the component contribution estimation process.
        Returns:
            the analysis results as a list of dictionaries

    is_using_group_contribution(reaction: equilibrator_api.phased_reaction.PhasedReaction) → bool
        Check whether group contribution is needed to get this reactions' dG.
        Returns:
            true iff group contribution is needed

    multicompartmental_standard_dg_prime(reaction_inner: equilibrator_api.phased_reaction.PhasedReaction, reaction_outer: equilibrator_api.phased_reaction.PhasedReaction, e_potential_difference: equilibrator_api.Q_, p_h_outer: equilibrator_api.Q_, ionic_strength_outer: equilibrator_api.Q_, p_mg_outer: equilibrator_api.Q_ = default_physiological_p_mg, tolerance: float = 0.0) → equilibrator_api.ureg.Measurement
        Calculate the transformed energies of a multi-compartmental reaction.
        Based on the equations from Harandsdottir et al. 2012 ([https://doi.org/10.1016/j.bpj.2012.02.032](https://doi.org/10.1016/j.bpj.2012.02.032))
        Parameters:
            reaction_inner (PhasedReaction) – the inner compartment half-reaction
            reaction_outer (PhasedReaction) – the outer compartment half-reaction
            e_potential_difference (Quantity) – the difference in electro-static potential between the outer and inner compartments
            p_h_outer (Quantity) – the pH in the outside compartment
            ionic_strength_outer (Quantity) – the ionic strength outside
            p_mg_outer (Quantity (optional)) – the pMg in the outside compartment
            tolerance (Float (optional)) – tolerance for identifying inbalance between inner and outer reactions (default = 0)
        Returns:
            standard_dg_prime – the transport reaction Gibbs free energy change

    static parse_formula_side(s: str) → Dict[str, float]
        Parse one side of the reaction formula.

    static parse_formula(formula: str) → Dict[str, float]
        Parse a two-sided reaction formula.

    static create_stoichiometric_matrix_from_reaction_formulas(formulas: Iterable[str]) → pandas.DataFrame
        Build a stoichiometric matrix.
        Parameters:
            formulas (Iterable[str]) – String representations of the reactions.
        Returns:
            The stoichiometric matrix as a DataFrame whose indexes are the compound IDs and its columns are the reaction IDs (in the same order as the input).

    create_stoichiometric_matrix_from_reaction_objects(reactions: Iterable[equilibrator_api.phased_reaction.PhasedReaction]) → pandas.DataFrame
        Build a stoichiometric matrix.
        Parameters:
            reactions (Iterable[PhasedReaction]) – The collection of reactions to build a stoichiometric matrix from.
        Returns:
            The stoichiometric matrix as a DataFrame whose indexes are the compounds and its columns are the reactions (in the same order as the input).
```

# equilibrator_api.phased_compound

inherit from equilibrator_cache.models.compound.Compound an add phases.

**Module Contents**
```python
equilibrator_api.phased_compound.AQUEOUS_PHASE_NAME = aqueous

equilibrator_api.phased_compound.GAS_PHASE_NAME = gas

equilibrator_api.phased_compound.LIQUID_PHASE_NAME = liquid

equilibrator_api.phased_compound.SOLID_PHASE_NAME = solid

equilibrator_api.phased_compound.REDOX_PHASE_NAME = redox

equilibrator_api.phased_compound.PhaseInfo

equilibrator_api.phased_compound.PHASE_INFO_DICT

equilibrator_api.phased_compound.NON_AQUEOUS_COMPOUND_DICT

equilibrator_api.phased_compound.MicroSpecie

equilibrator_api.phased_compound.PHASED_COMPOUND_DICT

equilibrator_api.phased_compound.CARBONATE_INCHIS
```
```python
class equilibrator_api.phased_compound.Condition(phase: str, abundance: equilibrator_api.ureg.Quantity = None)
    Bases: object
    A class for defining the conditions of a compound.
    I.e. the phase and the abundance.

    property phase → str
        Return the phase.

    property abundance → equilibrator_api.ureg.Quantity
        Return the abundance.

    property standard_abundance → equilibrator_api.ureg.Quantity
        Return the standard abundance in this phase.

    property physiolical_abundance → equilibrator_api.ureg.Quantity
        Return the default physiological abundance in this phase.

    property dimensionality → str
        Return the dimensionality of the abundance in this phase.
        E.g. [concentration] for aqueous phase, or [pressure] for gas phase.
        Returns:
            the dimensionality in this phase, or None if abundance is fixed.

    property ln_abundance → float
        Return the log of the ratio between given and std abundances.

    property ln_physiological_abundance → float
        Return the log of the ratio between phys and std abundances.

    reset_abundance() → None
        Reset the abundance to standard abundance.

    property is_physiological → bool
        Return True iff the abundance is the same as the physiological.
        Returns:
            True if the abundance is in physiological conditions, or if the abundance if fixed in this phase anyway.
```
```python
class equilibrator_api.phased_compound.PhasedCompound(compound: equilibrator_cache.Compound, condition: Condition = None)
    Bases: object
    A class that combines a equilibrator_api Compound and a Condition.

    static get_default(compound: equilibrator_cache.Compound) → Condition
        Get the default phase of a compound.
        Parameters:
            compound – a Compound
        Returns:
            the default phase

    property atom_bag → dict
        Get the compound's atom bag.

    property smiles → str
        Get the compound's InChI.

    property inchi → str
        Get the compound's InChI.

    property inchi_key → str
        Get the compound's InChIKey.

    property id → int
        Get the compound's equilibrator internal ID.

    property formula → str
        Get the chemical formula.

    property mass → float
        Get the chemical molecular mass.

    property phase → str
        Get the phase.

    property html_formula → str
        Get the chemical formula.

    property phase_shorthand → str
        Get the phase shorthand (i.e. 'l' for liquid).

    property possible_phases → Tuple[str]
        Get the possible phases for this compound.

    property abundance → equilibrator_api.ureg.Quantity
        Get the abundance.

    property ln_abundance → float
        Return the log of the abundance (for thermodynamic calculations).

    property ln_physiological_abundance → float
        Return the log of the default physiological abundance.

    property is_physiological → bool
        Check if the abundance is physiological.

    get_stored_standard_dgf_prime(p_h: equilibrator_api.ureg.Quantity, ionic_strength: equilibrator_api.ureg.Quantity, temperature: equilibrator_api.ureg.Quantity, p_mg: equilibrator_api.ureg.Quantity) → equilibrator_api.ureg.Quantity
        Return the stored formation energy of this phased compound.
        Only if it exists, otherwise return None (and we will use component-contribution later to get the reaction energy).
        Parameters:
            p_h –
            ionic_strength –
            temperature –
            p_mg –
        Returns:
            standard_dgf_prime (in kJ/mol)

    get_stored_standard_dgf() → equilibrator_api.ureg.Quantity
        Return the stored formation energy of this phased compound.
        Only if it exists, otherwise return None (and we will use component-contribution later to get the reaction energy).
        Returns:
            standard_dgf (in kJ/mol)

    get_stored_microspecie() → MicroSpecie
        Get the stored microspecies (from the PHASED_COMPOUND_DICT).
        Returns:
            The MicroSpecie namedtuple with the stored formation energy, or None if this compound has no stored value at this phase.

    serialize() → dict
        Return a serialized version of all the compound thermo data.
        Returns:
            a list of dictionaries with all the microspecies data
```
```python
class equilibrator_api.phased_compound.Proton(compound: equilibrator_cache.Compound)
    Bases: PhasedCompound
    A class specifically for protons.

    property abundance → equilibrator_api.ureg.Quantity
        Get the abundance.

    property ln_physiological_abundance → float
        Return the log of the default physiological abundance.

    property ln_abundance → float
        Return the log of the abundance (for thermodynamic calculations).
```
```python
class equilibrator_api.phased_compound.RedoxCarrier(compound: equilibrator_cache.Compound, potential: Optional[equilibrator_api.ureg.Quantity] = None)
    Bases: PhasedCompound
    A class specifically for redox carriers (with a given potential).

    get_stored_standard_dgf_prime(p_h: equilibrator_api.ureg.Quantity, ionic_strength: equilibrator_api.ureg.Quantity, temperature: equilibrator_api.ureg.Quantity, p_mg: equilibrator_api.ureg.Quantity) → equilibrator_api.ureg.Quantity
        Get the standard formation G'.

    get_stored_standard_dgf() → equilibrator_api.ureg.Quantity
        Get the standard formation G.

    property atom_bag → dict
        Get the compound's atom bag.

    property ln_abundance → float
        Return the log of the abundance (for thermodynamic calculations).

    property ln_physiological_abundance → float
        Return the log of the default physiological abundance.

    property is_physiological → bool
        Check if the abundance is physiological.
```

# equilibrator_api.phased_reaction

inherit from equilibrator_cache.reaction.Reaction an add phases.

**Module Contents**
```python
class equilibrator_api.phased_reaction.PhasedReaction(sparse: Dict[equilibrator_cache.Compound, float], arrow: str = '<=>', rid: str = None, sparse_with_phases: Dict[equilibrator_api.phased_compound.PhasedCompound, float] = None)
    Bases: equilibrator_cache.Reaction
    A daughter class of Reaction that adds phases and abundances.

    REACTION_COUNTER = 0

    static to_phased_compound(cpd: equilibrator_cache.Compound) → equilibrator_api.phased_compound.PhasedCompound
        Convert a Compound object to PhasedCompound.

    clone() → PhasedReaction
        Clone this reaction object.

    reverse() → PhasedReaction
        Return a PhasedReaction with the reverse reaction.

    get_element_data_frame() → pandas.DataFrame
        Tabulate the elemental composition of all reactants.
        Returns:
            A data frame where the columns are the compounds and the indexes are atomic elements.

    hash_md5(reversible: bool = True) → str
        Return a MD5 hash of the PhasedReaction.
        This hash is useful for finding reactions with the exact same stoichiometry. We create a unique formula string based on the Compound IDs and coefficients.
        Parameters:
            reversible (bool) – a flag indicating whether the directionality of the reaction matters or not. if True, the same value will be returned for both the forward and backward versions.
        Returns:
            hash – a unique hash string representing the Reaction.

    set_abundance(compound: equilibrator_cache.Compound, abundance: equilibrator_api.ureg.Quantity)
        Set the abundance of the compound.

    reset_abundances()
        Reset the abundance to standard levels.

    set_phase(compound: equilibrator_cache.Compound, phase: str)
        Set the phase of the compound.

    get_phased_compound(compound: equilibrator_cache.Compound) → Tuple[equilibrator_api.phased_compound.PhasedCompound, float]
        Get the PhasedCompound object by the Compound object.

    get_phase(compound: equilibrator_cache.Compound) → str
        Get the phase of the compound.

    get_abundance(compound: equilibrator_cache.Compound) → equilibrator_api.ureg.Quantity
        Get the abundance of the compound.

    property is_physiological → bool
        Check if all concentrations are physiological.
        This function is used by eQuilibrator to know if to present the adjusted dG' or not (since the physiological dG' is always displayed and it would be redundant).
        Returns:
            True if all compounds are at physiological abundances.

    get_stoichiometry(compound: equilibrator_cache.Compound) → float
        Get the abundance of the compound.

    add_stoichiometry(compound: equilibrator_cache.Compound, coeff: float) → None
        Add to the stoichiometric coefficient of a compound.
        If this compound is not already in the reaction, add it.

    separate_stored_dg_prime(p_h: equilibrator_api.ureg.Quantity, ionic_strength: equilibrator_api.ureg.Quantity, temperature: equilibrator_api.ureg.Quantity, p_mg: equilibrator_api.ureg.Quantity) → Tuple[equilibrator_cache.Reaction, equilibrator_api.ureg.Quantity]
        Split the PhasedReaction to aqueous phase and all the rest.
        Parameters:
            p_h –
            ionic_strength –
            temperature –
            p_mg –
        Returns:
            a tuple (residual_reaction, stored_dg_prime) where residual_reaction is a Reaction object (excluding the compounds that had stored values), and stored_dg_prime is the total dG' of the compounds with stored values (in kJ/mol).

    separate_stored_dg() → Tuple[equilibrator_cache.Reaction, equilibrator_api.ureg.Quantity]
        Split the PhasedReaction to aqueous phase and all the rest.
        Returns:
            a tuple (residual_reaction, stored_dg) where residual_reaction is a Reaction object (excluding the compounds that had stored values), and stored_dg is the total dG of the compounds with stored values (in kJ/mol).

    dg_correction() → equilibrator_api.ureg.Quantity
        Calculate the concentration adjustment in the dG' of reaction.
        Returns:
            the correction for delta G in units of RT

    physiological_dg_correction() → equilibrator_api.ureg.Quantity
        Calculate the concentration adjustment in the dG' of reaction.
        Assuming all reactants are in the default physiological concentrations (i.e. 1 mM)
        Returns:
            the correction for delta G in units of RT

    serialize() → List[dict]
        Return a serialized version of all the reaction thermo data.
```

# equilibrator_api.reaction_parser

A parser for reaction formulae.

**Module Contents**
```python
equilibrator_api.reaction_parser.POSSIBLE_REACTION_ARROWS = ['<=>', '<->', '-->', '<--', '=>', '<=', '->', '<-', '=', '⇌', '→', '←', '⇒']

equilibrator_api.reaction_parser.make_reaction_parser() → pyparsing.Forward
    Build pyparsing-based recursive descent parser for chemical reactions.
    Returns:
        parser
```
```python
equilibrator_api.default_phase = aqueous

equilibrator_api.default_physiological_p_h

equilibrator_api.default_physiological_p_mg

equilibrator_api.default_physiological_ionic_strength

equilibrator_api.default_physiological_temperature

equilibrator_api.default_conc_lb

equilibrator_api.default_conc_ub

equilibrator_api.default_e_potential

equilibrator_api.default_rmse_inf
```