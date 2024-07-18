import numpy as npfrom scipy.optimize import curve_fitfrom scipy.signal import convolvefrom ..core import TACfrom .kineticmodel import ReferenceKineticModeldef model(reftac: TAC, R1: float, k2: float, BPND: float):        """    reftac : reference tissue tac    R1 : K1/K1p    BPND : k3/k4    """        k2a = k2/(1+BPND)            t = reftac.frameschedule.mid_points        t_upsampled, dt = np.linspace(np.min(t), np.max(t), 2048, retstep=True)    # print(f't_upsampled: {t_upsampled}')    # print(f'dt: {dt}')        if reftac.data.shape[0] != 1:        raise ValueError("The reftac.data.shape[0] should be 1")        reftac_upsampled = np.interp(t_upsampled, t, reftac.data.flatten())        conv_upsampled = convolve(reftac_upsampled, np.exp(-k2a * t_upsampled), mode='full')[:len(t_upsampled)] * dt        tac_upsampled = R1 * reftac_upsampled + (k2 - R1 * k2a) * conv_upsampled        tac = np.interp(t, t_upsampled, tac_upsampled)    return tacclass SRTM_Model(ReferenceKineticModel):    def __init__(self,                  reftac: TAC,                  tacs: TAC):                super().__init__(reftac, tacs)                self.micro_params = {'R1': None,                            'k2': None,                            'BPND': None}        self.macro_params = {}        self.param_unit = {'R1': 'unitless',                           'k2': '/min',                           'BPND': 'unitless'}    def fit(self):                R1_arr = np.zeros(self.tacs.num_elements)        k2_arr = np.zeros(self.tacs.num_elements)        BPND_arr = np.zeros(self.tacs.num_elements)                for i in range(self.tacs.num_elements):                                    pars, _ = curve_fit(model, self.reftac, self.tacs.data[i,:])                    R1, k2, BPND = pars                        R1_arr[i] = R1            k2_arr[i] = k2            BPND_arr[i] = BPND                    self.set_parameter('R1', R1_arr, 'micro')        self.set_parameter('k2', k2_arr, 'micro')        self.set_parameter('BPND', BPND_arr, 'micro')                return None            def generate_fitted_tacs(self):                fitted_tacs_data = np.zeros((self.tacs.num_elements, self.tacs.frameschedule.num_frames))                            for i in range(self.tacs.num_elements):                    R1 = self.get_parameter('R1')[i]            k2 = self.get_parameter('k2')[i]            BPND = self.get_parameter('BPND')[i]                                        fitted_tacs_data[i,:] = model(self.reftac, R1, k2, BPND)                    self.fitted_tacs = TAC(frameschedule = self.tacs.frameschedule,                               data = fitted_tacs_data,                               rois = self.tacs.rois,                               unit = self.tacs.unit,)                        return None