import mibian
from scipy import stats
from numpy import log, exp, sqrt


class BSMcalc:

    """
        Calculates the price for a call or a put option

        Parameters:
            S (int or float): Stock Price
            K (int or float): Strike Price
            T (int or float): Time to expiration, measure in terms of years
            sigma (float): Standard Deviation
            rfr (float, optional): Risk-free Rate.
            D (float, optional): Annual Dividend yield.
            option_type (str, optional): Type of the option.
                                         Only accepts 'call' or 'put'.
            callprice (int or float): Price of the call option in order to
                                      calculate the implied volatility.
            putprice (int or float): Price of the put option in order to
                                     calculate the implied volatility.

        Returns:
            Price (float): Price of the option.
            nd1 (float): Probability of d1.
            nd2 (float): Probability of d2.
            implied_vol (float): Given the price of the option, it calculates
                                 the implied volatility in the option
    """

    def __init__(self, S, K, T, sigma, rfr=0.02, D=0.0, option_type='call'):
        self.S = S * exp(-D*T)
        self.K = K
        self.T = T
        self.sigma = sigma
        self.rfr = rfr
        self.D = D
        self.option_type = option_type
        self.d1 = ((log(self.S/self.K)+(self.rfr+self.sigma**2/2.)*self.T)
                   / (self.sigma*sqrt(self.T)))
        self.d2 = self.d1 - self.sigma*sqrt(self.T)

    def __str__(self):
        self.c = 'Stock price:{S} \nStrike price:{K} \nTime to maturity:{T} \
            \nVolatility:{sigma} \nRisk-free Rate:{rfr} \nDividends:{D} \
                \nOption type={op}'.format(S=self.S, K=self.K, T=self.T,
                                           sigma=self.sigma, rfr=self.rfr,
                                           D=self.D, op=self.option_type)
        return self.c

    def __repr__(self):
        self.c = 'Stock price:{S} \nStrike price:{K} \nTime to maturity:{T} \
            \nVolatility:{sigma} \nRisk-free Rate:{rfr} \nDividends:{D} \
                \nOption type={op}'.format(S=self.S, K=self.K, T=self.T,
                                           sigma=self.sigma, rfr=self.rfr,
                                           D=self.D, op=self.option_type)
        return self.c

    def price(self):

        """
            Calculates the price of the option.
        """

        if self.option_type == 'call':
            self.price = (self.S*stats.norm.cdf(self.d1)
                          - self.K*exp(self.rfr*-1*self.T)
                          * stats.norm.cdf(self.d2))
        elif self.option_type == 'put':
            self.price = (self.K*exp(-self.rfr*self.T)
                          * (1-stats.norm.cdf(self.d2))
                          - self.S*(1-stats.norm.cdf(self.d1)))
        else:
            raise ValueError('Only accepts call or put')
        return self.price

    def nd1(self):

        """
            Calculates how deep at-the-money, the option will expire.
        """

        if self.option_type == 'call':
            self.nd1 = stats.norm.cdf(self.d1)
        elif self.option_type == 'put':
            self.nd1 = stats.norm.cdf(self.d1) - 1
        else:
            raise ValueError('Input only call or put')
        return self.nd1

    def nd2(self):

        """
            Calculates the probability of the option to be exercised.
        """

        if self.option_type == 'call':
            self.nd2 = stats.norm.cdf(self.d2)
        elif self.option_type == 'put':
            self.nd2 = stats.norm.cdf(self.d2) - 1
        else:
            raise ValueError('Input only call or put')
        return self.nd2

    @staticmethod
    def implied_vol(S, K, T, rfr, callprice=0.0, putprice=0.0):

        """
            Calculates the implied volatility of the option.
        """

        S1 = float(S)
        K1 = float(K)
        interest_rate = float(rfr * 100)
        days_to_expiration = float(T * 365)
        if callprice > 0.0:
            bsc = mibian.BS([S1, K1, interest_rate,
                            days_to_expiration],
                            callPrice=callprice)
        elif putprice > 0.0:
            bsc = mibian.BS([S1, K1, interest_rate,
                            days_to_expiration],
                            putPrice=putprice)
        implied_vol = round(bsc.impliedVolatility/100, 6)
        return implied_vol


class Greeks(BSMcalc):

    """
        Calculates different types of sensitivities of the options.

        Parameters:
            S (int or float): Stock Price
            K (int or float): Strike Price
            T (int or float): Time to expiration, measure in terms of years
            sigma (float): Standard Deviation
            rfr (float, optional): Risk-free Rate.
            D (float, optional): Annual Dividend yield.
            option_type (str): Type of the option.
                                Only accepts 'call' or 'put'.

        Returns:
            delta (float): Rate of change of the option price with respect to
                            the price of the underlying asset.
            theta (float): Rate of change of the value of the portfolio
                            with respect to the passage of the time.
            gamma (float): Rate of change of the portfolio's delta
                            with respect to the price of the underlying asset.
            vega (float): Rate of change of the value of the portfolio
                       with respect to the volatility of the underlying asset.
            rho (float): Rate of change of the value of the portfolio
                            with respect to the interest rate.
    """

    def __init__(self, S, K, T, sigma, rfr=0.02, D=0.00, option_type='call'):
        BSMcalc.__init__(self, S, K, T, sigma, rfr=0.02, D=0.00,
                         option_type='call')
        self.price = BSMcalc.price(self)
        self.nd1 = BSMcalc.nd1(self)
        self.nd2 = BSMcalc.nd2(self)
        self.dnd1 = stats.norm.cdf((1/sqrt(2*3.14))*exp(-(self.d1*self.d1)/2))

    def __str__(self):
        return BSMcalc.__str__(self)

    def __repr__(self):
        return BSMcalc.__repr__(self)

    def delta(self):
        self.delta = self.nd1
        return self.delta

    def theta(self):
        if self.option_type == 'call':
            self.theta = (-(self.S*self.dnd1*self.sigma)/(2*sqrt(self.T))
                          - self.rfr*self.K*exp(-self.rfr*self.T)*self.nd2)
        elif self.option_type == 'put':
            self.theta = (-(self.S*self.dnd1*self.sigma)/(2*sqrt(self.T))
                          + self.rfr*self.K*exp(-self.rfr*self.T)*self.nd2)
        else:
            raise ValueError
        return self.theta

    def gamma(self):
        self.gamma = self.dnd1/(self.S*self.sigma*sqrt(self.T))
        return self.gamma

    def vega(self):
        self.vega = self.S*sqrt(self.T)*self.dnd1
        return self.vega

    def rho(self):
        if self.option_type == 'call':
            self.rho = self.K*self.T*exp(-self.rfr*self.T)*self.nd2
        elif self.option_type == 'put':
            self.rho = -self.K*self.T*exp(-self.rfr*self.T)*(self.nd2 - 1)
        else:
            raise ValueError
        return self.rho


class BrownRobinson(BSMcalc):

    """
        Calculates the option price considering the skewness and kurtosis of
        the underlaying asset.

        Parameters:
            S (int or float): Stock Price
            K (int or float): Strike Price
            T (int or float): Time to expiration, measure in terms of years
            sigma (float): Standard Deviation
            rfr (float, optional): Risk-free Rate.
            D (float, optional): Annual Dividend yield.
            option_type (str, optional): Type of the option.
                                        Only accepts 'call' or 'put'
            skewness (int or float, optional):
                                        Skewness of the underlying asset.
            kurtosis (int or float, optional):
                                        Kurtosis of the underlying asset.

        Returns:
            br_price (float): Price of the option.
    """

    def __init__(self, S, K, T, sigma, rfr=0.02, D=0.0, option_type='call',
                 skewness=0.0, kurtosis=3):
        BSMcalc.__init__(self, S, K, T, sigma, rfr=0.02, D=0.0,
                         option_type='call')
        self.cbs = BSMcalc.price(self)
        self.skewness = skewness
        self.kurtosis = kurtosis
        self.d = ((log(self.S/self.K)+(self.rfr+(self.sigma**2)/2)*self.T)
                  / (self.sigma*sqrt(self.T)))
        self.q3 = (0.1666*self.S*self.sigma*sqrt(self.T)
                   * ((2*self.sigma*sqrt(self.T))*stats.norm.cdf(self.d)
                   + self.sigma**2*self.T*stats.norm.cdf(self.d)))
        self.q4 = (0.04166*self.S*self.sigma*sqrt(self.T)
                   * ((self.d**2-1-3*self.sigma*sqrt(self.T)
                       * (self.d-self.sigma*sqrt(self.T)))
                   * stats.norm.cdf(self.d)
                   + self.sigma**3*self.T**1.5*stats.norm.cdf(self.d)))

    def __str__(self):
        self.d = ('Stock price:{S} \nStrike price:{K} \nTime to maturity:{T} \
            \nVolatility:{sigma} \nRisk-free Rate:{rfr} \nDividends:{D} \
                \nOption type={op} \nSkewness:{skew} \nKurtosis:{kurt} '
                  .format(S=self.S, K=self.K, T=self.T,
                          sigma=self.sigma, rfr=self.rfr, D=self.D,
                          op=self.option_type, skew=self.skewness,
                          kurt=self.kurtosis))
        return self.d

    def __repr__(self):
        self.d = ('Stock price:{S} \nStrike price:{K} \nTime to maturity:{T} \
            \nVolatility:{sigma} \nRisk-free Rate:{rfr} \nDividends:{D} \
                    \nOption type={op} \nSkewness:{skew} \nKurtosis:{kurt} '
                  .format(S=self.S, K=self.K, T=self.T,
                          sigma=self.sigma, rfr=self.rfr, D=self.D,
                          op=self.option_type, skew=self.skewness,
                          kurt=self.kurtosis))
        return self.d

    def br_price(self):
        self.cgc = self.cbs+self.skewness*self.q3+(self.kurtosis-1)*self.q4
        return self.cgc
