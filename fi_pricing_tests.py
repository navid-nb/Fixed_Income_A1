

from xml.parsers.expat import model
import numpy as np
from src.fi_pricing import NelsonSiegelSvensson, VasicekModel, CIRModel, HullWhiteModel

def assert_approx_equal(result, expected, tol=1e-3, message=""):
    assert np.isclose(result, expected, rtol=tol), \
        f"{message} failed \n \t\t Expected {expected}, got {result:.6f}"
    print(f"{message} passed: \n \t\t result: {result:.6f}, expected: {expected:.6f}")


def test_nss():
    pass



def test_vasicek():
    # Test1: Based on Quiz 1, question 1:
    model = VasicekModel(kappa=0.1, theta=0.06, sigma=0.15)
    result = 1000* model.P(t=0.0, T=3.0, rt=0.05)

    expected = 929.845
    assert_approx_equal(result, expected, message="vasicek test1")

    # Test2: Based on exercise 3.4:
    model = VasicekModel(kappa=0.1, theta=0.05, sigma=0.02)
    result = model.zcb_option(
        t=0.0,          # valuation time
        T=5.0,          # bond maturity
        rt=0.04,        # current short rate
        T_expiry=1.0,   # option expiry
        K=0.8,          # strike price
        option_type="call"
    )
    expected = 0.05115
    assert_approx_equal(result, expected, message="vasicek test2")


def test_CIRModel():
    # Test 1: Based on exercise 3.5:
    model = CIRModel(kappa=0.2, theta=0.04, sigma=0.1)
    #  Compute the value of the Treasury bond with 3-year maturity, 5% coupon rate and semi-annual coupon payments.
    result = model.coupon_bond_price(
        t=0.0,          # valuation time
        rt=0.06,        # current short rate
        cash_flows=np.array([2.5, 2.5, 2.5, 2.5, 2.5, 102.5]), 
        payment_dates=np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0]),
    )
    expected = 98.5352
    assert_approx_equal(result, expected, message="CIRModel test1")

    # result = model.coupon_bond_option(
    #     t=0.0,          # valuation time
    #     rt=0.06,        # current short rate
    #     T_expiry=1.0,   # option expiry
    #     K_bond=99.5,    # strike price
    #     cash_flows=np.array([2.5, 2.5, 2.5, 2.5, 2.5, 102.5]), 
    #     payment_dates=np.array([0.5, 1.0, 1.5, 2.0, 2.5, 3.0]),
    #     option_type="call"
    # )
    # expected = 1.02683
    # # assert_approx_equal(result, expected, message="CIRModel test2")

    # Test 2: Based on quiz 1 question 2:
    model = CIRModel(kappa=0.1, theta=0.06, sigma=0.15)
    result = 1000* model.zcb_option(
        t=0.0,          # valuation time
        T=5.0,          # bond maturity
        rt=0.05,        # current short rate
        T_expiry=2.0,   # option expiry
        K=model.P(t=0, T=5.0, rt=0.05),          # strike price (ATM)
        option_type="call"
    )
    expected = 84.895 
    assert_approx_equal(result, expected, message="CIRModel test2")


def test_HullWhite():
    # Test 1 & 2: Based on Quiz 1 Question 3 & 4:
    # Create NSS curve
    nss = NelsonSiegelSvensson(a=0.035, b=0.015, c=0.8,
        d=-0.7, tau=1.5, theta=3.0)
    
    r0 = nss.zcy(t=0, T=0)
    model = HullWhiteModel(kappa=0.1, sigma=0.1, yield_curve=nss)
    cash_flows, payment_dates = model.coupon_bond_cashflow_calculator(coupon=0.05, maturity=10, frequency=2, face_value=1000)

    result1 = model.coupon_bond_price(
        t=3/365,
        cash_flows=cash_flows,
        payment_dates=payment_dates,
        rt=r0
    )
    
    expected1 = 1734.85
    assert_approx_equal(result1, expected1, message="HullWhite test1")

    result2 = model.coupon_bond_option(
        t=3/365,
        rt=r0,
        T_expiry=3,
        K_bond=1600,
        cash_flows=cash_flows,
        payment_dates=payment_dates,
        option_type="call"
    )
    expected2 = 614.278
    assert_approx_equal(result2, expected2, message="HullWhite test2")


if __name__ == "__main__":
    test_nss()
    test_vasicek()
    test_CIRModel()
    test_HullWhite()
