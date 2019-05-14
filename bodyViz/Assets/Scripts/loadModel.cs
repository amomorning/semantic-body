using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;
using System.IO;

public class loadModel : MonoBehaviour
{
    void Start()
    {
        Button btn = this.GetComponent<Button>();
        btn.onClick.AddListener(OnClick);
    }

    private void OnClick()
    {
        var del = GameObject.Find("TEMP(Clone)");
        if (del != null)
        {
            Destroy(del);
        }
        var ave = GameObject.Find("AVE(Clone)");
        if (ave != null) GameObject.Destroy(ave);

        var prefab = Resources.Load("TEMP");
        var go = UnityEngine.Object.Instantiate(prefab, UnityEngine.Vector3.zero, Quaternion.identity) as GameObject;
        if (go == null)
        {
            Debug.Log("why ??????????");
        }
        else
        {
            Debug.Log(go.name);
            var defaultGo = go.transform.FindChild("default").gameObject;
            var comp = defaultGo.AddComponent<RotateTemp>();
            comp.AVE = defaultGo;
            Debug.Log(comp.AVE == null ? "AVE NULL !!!!" : comp.AVE.name);
        }
        //GameObject ave = GameObject.Find("AVE");
        //if (ave == null)
        //{
        //    Debug.Log("AVE is null");
        //}
        //else ave.SetActive(false);
        //Debug.Log("hhhhhhhhhhhhhhhhh");
        ////GameObject tmp = GameObject.Find("TEMP").GetComponentInChildren<GameObject>();
        //GameObject tmp = GameObject.Find("TEMP");
        //var it = tmp.transform.FindChild("default");
        //if (it == null)
        //{
        //    Debug.Log("??????");
        //}
        //it.gameObject.SetActive(true);
        //tmp.SetActive(true);

        Debug.Log("Button Clicked. ClickHandler.");
    }

}
