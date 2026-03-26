'''
Code       : off_module.py
Date       : 2024.11.01
Contributer: H.Y.Li (liyuhua0909@126.com), G.Y.Chen (gychen@smail.nju.edu.cn)
Purpose    : Build a Outflow Field magnetic field model...

### --------------------------------- ###
Remark:
2024.11.01: Build the code
2024.11.22: Modified the rc calculation
'''

from .needs import *
from scipy.special import lambertw
import pkg_resources
import pickle
from scipy import linalg, interpolate

from .pfss_module import pfss_solver
from .funcs import brtp2bxyz, trilinear_interpolation, Brtp_lm, Associated_Legendre,DAssociated_Legendre
from .magline import rk45, magline_stepper, magline_solver, show_boundary, show_maglines,parallel_magline_solver,show_current_sheet

# ==========================================

L0  = 6.995e10
t0  = 5972.5794
v0  = 1.16448846777562e7
nu0 = 5.e-17*L0**2/t0
GM  = 9.54e4
shooting_file = pkg_resources.resource_filename('pfss', 'shooting_target.pkl')
with open(shooting_file, 'rb') as f:
    shooting_target = pickle.load(f)

def interp_initial(v,l=0):
    vlist = shooting_target['vlist']
    init  = shooting_target['target'][l]
    interp1 = interpolate.interp1d(
        vlist, init[:,0], kind='cubic',copy=False, assume_sorted=True
    )
    interp2 = interpolate.interp1d(
        vlist, np.log10(init[:,1]), kind='cubic',copy=False, assume_sorted=True
    )
    ret = np.array([interp1(v), np.power(10,interp2(v))])
    return ret

def rk45_solver(rfun, x, y, dl, **kwargs):
    sig = kwargs.get('sig', 1)
    x0  = x
    k1  = rfun(x0         , y        , **kwargs)
    k2  = rfun(x0+sig*dl/2, y+k1*dl/2, **kwargs)
    k3  = rfun(x0+sig*dl/2, y+k2*dl/2, **kwargs)
    k4  = rfun(x0+sig*dl  , y+k3*dl  , **kwargs)
    k   = (k1+2*k2+2*k3+k4)/6*sig
    ret = y+k*dl
    return ret

# ==========================================

class off_solver(pfss_solver):
    def __init__(self,
                 fits_file=None,
                 n_r      = 400,
                 n_t      = 200,
                 n_p      = 400,
                 lmax     = 80,
                 Rs       = 2.5,
                 v1       = 50
                ):
        super().__init__(fits_file,n_r,n_t,n_p,lmax,Rs)
        self.scs_file  = './Brtp_off.npy'
        self.save_name = 'off_solver.pkl'
        self.v1        = v1*1e5/v0
        self.rc        = self.cal_rc()
        self.off_file  = './Outflow_field.npy'

    def set_v1(self, v1):
        self.v1 = v1/v0*1e5

    # def cal_rc(self):
    #     r1 = self.Rs
    #     v1 = self.v1*v0*1e-5
    #     kk = 9.54e4
    #     e  = np.exp(1)
    #     rc = -3/4*r1*lambertw(-4/3*(r1*v1**2/e**4/kk)**(1/3),k=-1)
    #     return np.real(rc)

    # def vout(self, rho):
    #     v1  = self.v1
    #     r1  = self.Rs
    #     rc  = self.rc
    #     ret = np.exp(-2*(np.exp(-rho)*rc-rc/r1+rho))*r1**2*v1
    #     return ret
    
    # def d1_vout(self, rho):
    #     v1  = self.v1
    #     r1  = self.Rs
    #     rc  = self.rc
    #     ret = np.exp(-2*(np.exp(-rho)*rc-rc/r1+rho))*r1**2*v1
    #     return ret

    def cal_rc(self):
        r1 = self.Rs
        v1 = self.v1*v0*1e-5
        r0 = GM/(2*v1**2)
        e  = np.exp(1)
        rc = np.where(r0>r1,
                      3/(r1/r0-4)*lambertw(1/3*(r1/r0-4)*(r1/r0)**(1/3)*np.exp(-1),k=-1)*r1,
                      3/(r1/r0-4)*lambertw(1/3*(r1/r0-4)*(r1/r0)**(1/3)*np.exp(-1),k= 0)*r1
                     )
        # rc = 3/(r1/r0-4)*lambertw(1/3*(r1/r0-4)*(r1/r0)**(1/3)*np.exp(-1),k=-1)*r1
        return np.real(rc)

    def D(self, r):
        rc = self.rc
        if isinstance(rc,np.ndarray):
            if rc.ndim==1:
                rc=rc[:,None]
        Dr = (rc/r)**4*np.exp(4*(1-rc/r)-1)
        return Dr

    def vout(self, rho):
        r   = np.exp(rho)
        rc  = self.rc
        cs  = np.sqrt(GM/2/self.rc)*1e5/v0
        if isinstance(cs,np.ndarray) and len(cs)>1:
            if cs.ndim==1:
                cs=cs[:,None]
                rc=rc[:,None]
        Dr  = self.D(r)
        # ret = cs*np.real(np.sqrt(-lambertw(-Dr, k=0))) if r<=rc else cs*np.real(np.sqrt(-lambertw(-Dr,k=-1)))
        ret = np.where(r<=rc, cs*np.real(np.sqrt(-lambertw(-Dr, k=0))), cs*np.real(np.sqrt(-lambertw(-Dr,k=-1))))
        return ret

    def d1_vout(self, rho):
        r   = np.exp(rho)
        rc  = self.rc
        Dr  = self.D(r)
        # w   = np.real(lambertw(-Dr, k=0)) if r<=rc else np.real(lambertw(-Dr, k=-1))
        w   = np.where(r<=rc, np.real(lambertw(-Dr, k=0)), np.real(lambertw(-Dr, k=-1)))
        cs  = np.sqrt(GM/2/self.rc)*1e5/v0
        dDr = -4*(1-rc/r)*Dr
        ret = cs*0.5*(-w)**0.5/Dr*w/(1-w)*dDr
        return ret
        

    def d2_vout(self, rho):
        v1  = self.v1
        r1  = self.Rs
        rc  = self.rc
        ret = 2*np.exp(-2*np.exp(-rho)*rc+2*rc/r1-4*rho)*r1**2*(2*np.exp(2*rho)-5*np.exp(rho)*rc+2*rc**2)*v1
        return ret

    # def rfun(self, rho, y=None, l=0, **kwargs):
    #     y0,y1 = y
    #     k1 = (3-nu0*np.exp(rho)*self.vout(rho))
    #     k2 = -(l*(l+1)-2+3*nu0*np.exp(rho)*self.vout(rho)+nu0*np.exp(rho)*self.d2_vout(rho))
    #     rf1 = -k1*y1-k2*y0
    #     rf0 = y1
    #     ret = np.array([rf0,rf1])
    #     return ret

    # def rfun(self, rho, y=None, l=0, **kwargs):
    #     y0,y1 = y
    #     k1  = (4-(2-nu0)*np.exp(rho)*self.vout(rho))
    #     k2  = -(l*(l+1)-3+(5-nu0)*np.exp(rho)*self.vout(rho)-(1-nu0)*np.exp(2*rho)*self.vout(rho)**2+np.exp(rho)*self.d1_vout(rho))
    #     rf1 = -k1*y1-k2*y0
    #     rf0 = y1
    #     ret = np.array([rf0,rf1])
    #     return ret

    def rfun(self, rho, y=None, l=0, **kwargs):
        y0,y1 = y
        k1 = (3-nu0*np.exp(rho)*self.vout(rho))
        k2 = -(l*(l+1)-2+3*nu0*np.exp(rho)*self.vout(rho)+nu0*np.exp(rho)*self.d1_vout(rho))
        rf1 = -k1*y1-k2*y0
        rf0 = y1
        ret = np.array([rf0,rf1])
        return ret

    def time_integration(self, initial, **kwargs):
        r1    = kwargs.get('r1',self.Rs)
        dl    = kwargs.get('dl',1e-3)
        l     = kwargs.get('l' , 0)
        ns    = kwargs.get('max_steps', 100000)
        rmax  = kwargs.get('rmax', None)
        rhoc  = np.log(r1)
        rho0  = rhoc
        rlist = [rho0]
        sol   = [initial]
        iter  = 0
        while rho0>0 and iter<=ns:
            x = rho0
            y = sol[-1]
            sol.append(rk45_solver(self.rfun,x,y,dl, sig=-1,l=l))
            rho0 = rho0-dl
            if rho0<0:
                rho0 = 0
                dl   = x-0
                sol[-1] = rk45_solver(self.rfun,0,y,dl, sig=-1, l=l)
            rlist.append(rho0)
        if rmax is None or rmax<=r1:
            rlist = np.array(rlist)
            sol   = np.array(sol)
            return rlist, sol
        else:
            sol_out = [initial]
            iter=0
            rho0 = rhoc
            rlist_out = [rho0]
            while rho0<np.log(rmax) or iter<=ns:
                x = rho0
                y = sol_out[-1]
                sol_out.append(rk45_solver(self.rfun, x, y, dl, sig=1, l=l))
                rho0 = rho0+dl
                rlist_out.append(rho0)
            rlist_ret = np.array(rlist[::-1]+rlist_out[1:])
            sol_ret   = np.array(sol[::-1]  +sol_out[1:])
            return rlist_ret, sol_ret

    def shooting_method(self, l=0, **kwargs):
        aim_i   = kwargs.get('aim_i', -10)
        aim_f   = kwargs.get('aim_f', 100)
        aim_N   = kwargs.get('aim_N', 100)
        aim_try = np.linspace(aim_i,aim_f,aim_N)
        v0      = self.vout(0)
        trial_shooting = []
        for ss in aim_try:
            initial = np.array([0,10**(-ss)])
            _,sol   = self.time_integration(initial, l=l)
            # trial_shooting.append(sol[-1].sum())
            trial_shooting.append(sol[-1,1]+(1-nu0*v0)*sol[-1,0])
        root = brentq(lambda x: interp1d(aim_try,np.log10(np.abs(trial_shooting)),kind='cubic')(x), aim_i, aim_f)
        initial = np.array([0,10**(-root)])
        # rlist,sol = self.time_integration(initial,l=l)
        return initial

    def compute_coefficient(self, **kwargs):
        self.Br = kwargs.get('bound', self.Br)
        t0 = time.time()
        lmax   = kwargs.get('lmax', self.lmax)
        device = kwargs.get('device', 'cuda' if torch.cuda.is_available() else 'cpu')
        device = torch.device(device)
        Alm = {il: {} for il in range(lmax + 1)} # coefficient for sin(m phi)
        Blm = {il: {} for il in range(lmax + 1)} # coefficient for cos(m phi)
        rr,tt,pp = self.get_rtp()
        tt       = torch.from_numpy(tt[0]  ).to(device)
        pp       = torch.from_numpy(pp[0]  ).to(device)
        Br       = torch.from_numpy(self.Br).to(device)
        for m in range(0,lmax+1):
            for l in range(m,lmax+1):
                if l==m:
                    P_l00 = Associated_Legendre(l,  m, torch.cos(tt))
                    P_lp1 = Associated_Legendre(l+1,m, torch.cos(tt))
                    cosmp = torch.cos(m*pp)
                    sinmp = torch.sin(m*pp)
                    if l==0:
                        Alm[l][m] = 0
                        Blm[l][m] = (torch.sum(Br*P_l00*torch.sin(tt))/torch.sum(P_l00**2*torch.sin(tt))).item()
                    else:
                        Alm[l][m] = (torch.sum(Br*P_l00*sinmp*torch.sin(tt))/torch.sum(P_l00**2*sinmp**2*torch.sin(tt))).item()
                        Blm[l][m] = (torch.sum(Br*P_l00*cosmp*torch.sin(tt))/torch.sum(P_l00**2*cosmp**2*torch.sin(tt))).item()
                    P_ln1 = P_l00
                    P_l00 = P_lp1
                else:
                    P_lp1 = Associated_Legendre(l+1, m, torch.cos(tt), pn2=P_ln1, pn1=P_l00)
                    cosmp = torch.cos(m*pp)
                    sinmp = torch.sin(m*pp)
                    Alm[l][m] = (torch.sum(Br*P_l00*sinmp*torch.sin(tt))/torch.sum(P_l00**2*sinmp**2*torch.sin(tt))).item()
                    Blm[l][m] = (torch.sum(Br*P_l00*cosmp*torch.sin(tt))/torch.sum(P_l00**2*cosmp**2*torch.sin(tt))).item()
                    if m==0:
                        Alm[l][m]=0
                    P_ln1 = P_l00
                    P_l00 = P_lp1
        self.Alm = Alm
        self.Blm = Blm
        ti = time.time()
        print(f'Finishing computing coefficient takes {(ti-t0):8.3f} sec...')
        return Alm, Blm

    # def compute_Hl(self,r,l=0,**kwargs):
    #     r_list, sol = self.shooting_method(l=l,**kwargs)
    #     interp_func = interp1d(r_list, sol[:,0], kind='cubic')
    #     Hl          = interp_func(np.log(r))
    #     return Hl

    # def compute_Gl(self, r, l=0, **kwargs):
    #     Hl    = kwargs.get('Hl', self.compute_Hl(r,l=l))
    #     Gl    = np.zeros_like(Hl)
    #     Gl[0] = 1.
    #     for i in range(len(Hl)-1):
    #         Gl[i+1] = 0.5*l*(l+1)*Hl[i]*(1-r[i]**2/r[i+1]**2)+Gl[i]*r[i]**2/r[i+1]**2
    #     return Gl

    # def _single_task(self, rls, l=0):
    #     Hl = self.compute_Hl(rls, l=l)
    #     Gl = self.compute_Gl(rls, l=l, Hl=Hl)
    #     return l, Hl, Gl
    
    def compute_HG(self,r,l=0,**kwargs):
        method = kwargs.get('HG_method', 'shooting')
        if (self.v1*v0/1e5<1) or (self.v1*v0/1e5>200) or method=='shooting':
            initial = self.shooting_method(l=l,**kwargs)
        else:
            initial = interp_initial(self.v1*v0/1e5,l=l)
        r_list, sol = self.time_integration(initial,l=l, rmax=r.max())
        interp_Hl   = interp1d(r_list, sol[:,0], kind='cubic')
        Gsol        = sol[:,1]+(1-nu0*np.exp(r_list)*self.vout(r_list))*sol[:,0]
        interp_Gl   = interp1d(r_list, Gsol    , kind='cubic')
        Hl          = interp_Hl(np.log(r))
        Gl          = interp_Gl(np.log(r))
        return Hl,Gl

    def _single_task(self, rls, l=0):
        Hl,Gl = self.compute_HG(rls,l=l)
        return l, Hl, Gl

    def radial_function(self, rls, **kwargs):
        lmax = self.lmax
        n_cores = kwargs.get('n_cores', 4)
        ret = []
        with ProcessPoolExecutor(max_workers=n_cores) as executor:
            futures = []
            for l in range(lmax+1):
                futures.append(executor.submit(self._single_task, rls, l))
            for future in as_completed(futures):
                ret.append(list(future.result()))
        return ret

    def get_outflow_field(self, **kwargs):
        '''
        Compute the outflow field using the shooting method.
        
        Parameters:
        -----------
        lmax: int, optional
            The maximum order of the spherical harmonics.
        devices: list, optional
            The devices to use for the computation.
        python: str, optional
            The python command to use for the computation.
        fast_mode: bool, optional
            Whether to use the fast mode for the computation.
        n_cores: int, optional
            The number of cores to use for the computation.
        bound: np.ndarray, optional
            The boundary condition for the outflow field.
        off_file: str, optional
            The file name to save the outflow field.
        v_mode: str, optional
            The mode of the outflow field.
        v1: float, optional
            The initial velocity of the outflow field.
        rtp: tuple, optional
            The coordinates of the outflow field.
        skip: bool, optional
            Whether to skip the computation.

        Returns:
        --------
        np.ndarray
            The outflow field.
        '''
        t0       = time.time()
        device   = kwargs.get('device', 'cuda' if torch.cuda.is_available() else 'cpu')
        device   = torch.device(device)
        lmax     = kwargs.get('lmax', self.lmax)
        Alm      = kwargs.get('Alm', None)
        Blm      = kwargs.get('Blm', None)
        n_cores  = kwargs.get('n_cores', 50)
        FM       = kwargs.get('fast_mode', False)
        OF       = kwargs.get('off_file' , self.off_file)
        vM       = kwargs.get('v_mode', 'const') # 'const', 'theta', 'theta_phi'
        v1       = kwargs.get('v1' , None)
        if Alm is None or Blm is None:
            Alm,Blm = self.compute_coefficient(**kwargs)
        rtp      = kwargs.get('rtp', self.get_rtp())
        skip     = kwargs.get('skip', False)
        rr,tt,pp = rtp
        rls      = rr[:,0,0]
        tls      = tt[0,:,0]
        pls      = pp[0,0,:]
        Nr,Nt,Np = self.n_r,self.n_t,self.n_p
        rr,tt,pp = torch.from_numpy(np.stack([rr,tt,pp])).to(device)
        br = torch.zeros_like(rr, dtype=torch.float64).to(device)
        bt = torch.zeros_like(tt, dtype=torch.float64).to(device)
        bp = torch.zeros_like(pp, dtype=torch.float64).to(device)
        if not skip:
            if vM=='const':
                HL = {il: {} for il in range(lmax+1)}
                GL = {il: {} for il in range(lmax+1)}
                with ProcessPoolExecutor(max_workers=n_cores) as executor:
                    futures = []
                    for l in range(lmax+1):
                        futures.append(executor.submit(self._single_task, rls, l))
                    for future in as_completed(futures):
                        l,hl,gl = future.result()
                        HL[l]   = hl
                        GL[l]   = gl
                        # HL[l]   = hl[:,np.newaxis,np.newaxis].repeat(Nt,axis=1).repeat(Np,axis=2)
                        # GL[l]   = gl[:,np.newaxis,np.newaxis].repeat(Nt,axis=1).repeat(Np,axis=2)
            elif vM=='theta':
                # HL = {il:np.zeros((len(rls), len(tls), len(pls))) for il in range(lmax+1)}
                # GL = {il:np.zeros((len(rls), len(tls), len(pls))) for il in range(lmax+1)}
                HL = {il:np.zeros((len(rls), len(tls))) for il in range(lmax+1)}
                GL = {il:np.zeros((len(rls), len(tls))) for il in range(lmax+1)}
                if len(v1)!=self.n_t:
                    raise TypeError(f"Shape of `v1` should be consistant with `nth`: {self.n_t}")
                else:
                    v1_init = self.v1
                    unit_v  = 1e5/v0
                    for i,iv1 in enumerate(v1):
                        it0     = time.time()
                        self.v1 = iv1*unit_v
                        self.rc = self.cal_rc()
                        with ProcessPoolExecutor(max_workers=n_cores) as executor:
                            futures = []
                            for l in range(lmax+1):
                                futures.append(executor.submit(self._single_task, rls, l,**kwargs))
                            for future in as_completed(futures):
                                l,hl,gl = future.result()
                                # HL[l][:,i,:]=hl[:,np.newaxis].repeat(len(pls),axis=1)
                                # GL[l][:,i,:]=gl[:,np.newaxis].repeat(len(pls),axis=1)
                                HL[l][:,i]=hl
                                GL[l][:,i]=gl
                        print(f"Radial Process: {i+1:4d}/{len(v1)} | dt: {time.time()-it0:8.3f} sec | Wall Time: {(time.time()-t0)/60:8.3f} min")
                    self.v1 = v1_init
                    
            elif vM=='theta_phi':
                HL = {il:np.zeros(len(tls), len(tls), len(pls)) for il in range(lmax+1)}
                GL = {il:np.zeros(len(tls), len(tls), len(pls)) for il in range(lmax+1)}
                if v1.shape!=(len(tls), len(pls)):
                    raise TypeError(f"Shape of `v1`({v1.shape}) should be consistant with ({len(tls)}, {len(pls)})...")
                else:
                    v1_init = self.v1
                    unit_v  = 1e5/v0
                    for i,vt in enumerate(v1):
                        for j, vtp in enumerate(vt):
                            self.v1 = vtp*unit_v
                            self.rc = self.cal_rc()
                            temp    = self.radial_function(rls, **kwargs)
                            for l, hl, gl in temp:
                                HL[l][:,i,j] = hl
                                GL[l][:,i,j] = gl
                    self.v1 = v1_init
            else:
                raise ValueError("`v_mode` only can choose 'const', 'theta' or 'theta_phi'...")
            print(f'Complete calculating the radial function Hl and Gl, wall_time: {(time.time()-t0):8.3f} sec...')
        if not FM:
            for m in range(0, lmax+1):
                for l in range(m, lmax+1):
                    alm = Alm[l][m]
                    blm = Blm[l][m]
                    if l==m:
                        P_l00 = Associated_Legendre(l,  m, torch.cos(tt))
                        P_lp1 = Associated_Legendre(l+1,m, torch.cos(tt))
                        cosmp = torch.cos(m*pp)
                        sinmp = torch.sin(m*pp)
                        Hl    = HL[l]
                        Gl    = GL[l]
                        Hl    = torch.from_numpy(Hl).to(device)
                        Gl    = torch.from_numpy(Gl).to(device)
                        Qlm   = P_l00
                        DQlm  = DAssociated_Legendre(l,m,torch.cos(tt),P_l00=P_l00,P_lp1=P_lp1)
                        Pm1   = sinmp
                        Pm2   = cosmp
                        Dpm1  = m*cosmp
                        Dpm2  =-m*sinmp
                        if l==0:
                            br+=blm*Gl*Qlm*Pm2
                        else:
                            br+=Gl*Qlm*(alm*Pm1+blm*Pm2)
                            bt+=-Hl*DQlm*torch.sin(tt)*(alm*Pm1+blm*Pm2)
                            bp+=Hl/torch.sin(tt)*Qlm*(alm*Dpm1+blm*Dpm2)
                        P_ln1 = P_l00
                        P_l00 = P_lp1
                    else:
                        P_lp1 = Associated_Legendre(l+1, m, torch.cos(tt), pn2=P_ln1, pn1=P_l00)
                        cosmp = torch.cos(m*pp)
                        sinmp = torch.sin(m*pp)
                        Hl    = HL[l]
                        Gl    = GL[l]
                        Hl    = torch.from_numpy(Hl).to(device)
                        Gl    = torch.from_numpy(Gl).to(device)
                        Qlm   = P_l00
                        DQlm  = DAssociated_Legendre(l,m,torch.cos(tt),P_l00=P_l00,P_lp1=P_lp1)
                        Pm1   = sinmp
                        Pm2   = cosmp
                        Dpm1  = m*cosmp
                        Dpm2  =-m*sinmp
                        br+=Gl*Qlm*(alm*Pm1+blm*Pm2)
                        bt+=-Hl*DQlm*torch.sin(tt)*(alm*Pm1+blm*Pm2)
                        bp+=Hl/torch.sin(tt)*Qlm*(alm*Dpm1+blm*Dpm2)
                        P_ln1 = P_l00
                        P_l00 = P_lp1           
            br  =  np.real(br.detach().cpu().numpy())
            bt  =  np.real(bt.detach().cpu().numpy())
            bp  =  np.real(bp.detach().cpu().numpy())
            ret =  np.stack([br,bt,bp], axis=0)
        else:
            save_name = kwargs.get('save_name', self.save_name)
            devices   = kwargs.get('devices'  , ['cuda:0']    )
            PY        = kwargs.get('python'   , 'python '     )
            self.save_name  = save_name
            self.lmax       = lmax
            if not skip:
                self.info['HL'] = HL
                self.info['GL'] = GL
                self.info['rtp']= rtp
                self.info['v_mode']=vM
            # print('HL: ', HL[0].shape)
            nD = len(devices)
            di = (lmax+1)*(lmax+2)//2/nD
            m_assigned = [0]
            icnt = 0
            for m in range(lmax+1):
                icnt+=lmax+1-m
                if icnt > di:
                    m_assigned.append(m+1)
                    icnt = icnt % di
            m_assigned.append(lmax+1)
            self.info['mlist']   = m_assigned
            self.info['devices'] = devices
            self.save(save_name=save_name)
            commands = []
            for i,D in enumerate(devices):
                command = PY+f' -u -m pfss.off_scripts -i {save_name} -n {i}'
                commands.append(command)
            processes = []
            for command in commands:
                process = subprocess.Popen(command, shell=True)
                processes.append(process)
            for process in processes:
                process.wait()
            off_files = sorted(glob.glob('./OFF_temp_files/*.npy'))
            ret = np.zeros((3,Nr,Nt,Np))
            for f in off_files:
                ret += np.load(f)
            br,bt,bp = ret
            shutil.rmtree('./OFF_temp_files')
    
        self.Br_BB = br[ 0]
        self.Br_SS = br[-1]
        ti = time.time()
        self.off_file = OF
        np.save(OF, ret)
        print(f'Building the Outflow Field successfully using time: {(ti-t0)/60:8.3f} min...')
        return ret

    def save_vts(self, **kwargs):
        vts_name = kwargs.pop('vts_name', 'off')
        Brtp     = kwargs.pop('Brtp', np.load(self.off_file))
        super().save_vts(vts_name=vts_name, Brtp=Brtp, **kwargs)

    def save_vtu(self, **kwargs):
        vtu_name = kwargs.pop('vtu_name', 'off')
        Brtp     = kwargs.pop('Brtp', np.load(self.off_file))
        super().save_vtu(vtu_name=vtu_name, Brtp=Brtp,**kwargs)

    def show_magline(self, **kwargs):
        load_file = self.off_file
        super().show_magline(load_file=load_file, **kwargs)